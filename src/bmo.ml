open Owl

module D = struct
  type arr = Arr.arr

  let print_dim x =
    let shp = Arr.shape x in
    Array.iter (fun s -> Printf.printf "%i " s) shp;
    print_newline ();
    flush_all ()


  (* batch matrix multiplication *)
  let bmm =
    let check_dims_match ndim shpx shpy =
      if shpx.(ndim - 1) <> shpy.(ndim - 2)
      then failwith "bmm: last two dimensions mismatch"
      else
        for i = 0 to ndim - 3 do
          if shpx.(i) <> shpy.(i) then failwith "bmm: do not support broadcast"
        done
    in
    fun x y ->
      let shpx = Arr.shape x in
      let shpy = Arr.shape y in
      let ndimx = Array.length shpx in
      let ndimy = Array.length shpy in
      if ndimx <> ndimy then failwith "bmm: dimensions of [x] and [y] must be the same";
      if ndimx < 2
      then failwith "bmm: dimensions must be greater than 2"
      else if ndimx = 2
      then Arr.dot x y
      else (
        let ndim = ndimx in
        check_dims_match ndim shpx shpy;
        let shp = Array.copy shpx in
        let m = shpx.(ndim - 2) in
        let k = shpx.(ndim - 1) in
        let l = shpy.(ndim - 2) in
        assert (k == l);
        let n = shpy.(ndim - 1) in
        shp.(ndim - 1) <- n;
        shp.(ndim - 2) <- m;
        let batch_size = Array.fold_left ( * ) 1 shp / m / n in
        let z = Arr.empty shp in
        (* Printf.printf "%i, %i, %i %i\n%!" m n k batch_size; *)
        let xr = Arr.reshape x [| batch_size; m; k |] in
        let yr = Arr.reshape y [| batch_size; k; n |] in
        let zr = Arr.reshape z [| batch_size; m; n |] in
        for j = 0 to batch_size - 1 do
          let x1 = Bigarray.Genarray.sub_left xr j 1 in
          let x2 = Bigarray.Genarray.sub_left yr j 1 in
          let x3 = Bigarray.Genarray.sub_left zr j 1 in
          let alpha = 1. in
          let beta = 0. in
          let a = Arr.flatten x1 |> Bigarray.array1_of_genarray in
          let b = Arr.flatten x2 |> Bigarray.array1_of_genarray in
          let c = Arr.flatten x3 |> Bigarray.array1_of_genarray in
          let layout = Owl_cblas_basic.CblasRowMajor in
          let transa = Owl_cblas_basic.CblasNoTrans in
          let transb = Owl_cblas_basic.CblasNoTrans in
          Owl_cblas_basic.gemm layout transa transb m n k alpha a k b n beta c n
        done;
        z)


  let bchol ?(upper = true) x =
    let shp = Arr.shape x in
    let ndims = Array.length shp in
    if ndims < 2
    then failwith "bchol: dimension must be greater than 2"
    else if ndims = 2
    then Linalg.D.chol ~upper x
    else (
      let x = Arr.copy x in
      let m = shp.(ndims - 1) in
      let n = shp.(ndims - 2) in
      if m <> n
      then failwith "bchol: last two dimensions do not match"
      else (
        let batch_size = Array.fold_left ( * ) 1 shp / m / n in
        let xr = Arr.reshape x [| batch_size; m; n |] in
        for j = 0 to batch_size - 1 do
          let a =
            Bigarray.Genarray.sub_left xr j 1 |> fun x -> Arr.reshape x [| m; m |]
          in
          if upper
          then (
            Owl_lapacke.potrf ~uplo:'U' ~a |> ignore;
            for t = 1 to m - 1 do
              for s = 0 to t - 1 do
                Mat.set a t s 0.
              done
            done)
          else (
            Owl_lapacke.potrf ~uplo:'L' ~a |> ignore;
            for t = 1 to m - 1 do
              for s = 0 to t - 1 do
                Mat.set a s t 0.
              done
            done)
        done;
        x))
end

module AD = struct
  type t = Algodiff.D.t

  let rec _bmm =
    lazy
      (Algodiff.D.Builder.build_piso
         (module struct
           let label = "bmm"
           let ff_aa _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
           let ff_ab _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
           let ff_ba _ _ = raise Owl_exception.(NOT_IMPLEMENTED "bmm")
           let ff_bb a b = Algodiff.D.pack_arr D.(bmm a b)
           let df_da _cp _ap at bp = bmm at bp
           let df_db _cp ap _bp bt = bmm ap bt
           let df_dab _cp ap at bp bt = Algodiff.D.Maths.(bmm ap bt + bmm at bp)

           let dr_ab a b _cp ca =
             let open Algodiff.D in
             let dim_a = Array.length (shape a) in
             let dim_b = Array.length (shape b) in
             ( bmm !ca (Maths.swap (dim_b - 1) (dim_b - 2) (primal b))
             , bmm (Maths.swap (dim_a - 1) (dim_a - 2) (primal a)) !ca )


           let dr_a _a b _cp ca =
             let open Algodiff.D in
             let dim = Array.length (shape b) in
             bmm !ca (Maths.swap (dim - 1) (dim - 2) (primal b))


           let dr_b a _b _cp ca =
             let open Algodiff.D in
             let dim = Array.length (shape a) in
             bmm (Maths.swap (dim - 1) (dim - 2) (primal a)) !ca
         end : Algodiff.D.Builder.Piso))


  and bmm x = Stdlib.Lazy.force _bmm x
end

let%test _ =
  let module FD = Owl_algodiff_check.Make (Algodiff.D) in
  let n_samples = 1000 in
  let ff x =
    let x1 =
      Algodiff.D.Maths.get_slice [ []; [ 0; 1000 - 1 ] ] x
      |> fun x -> Algodiff.D.Maths.reshape x [| 10; 10; 10 |]
    in
    let x2 =
      Algodiff.D.Maths.get_slice [ []; [ 1000; -1 ] ] x
      |> fun x -> Algodiff.D.Maths.reshape x [| 10; 10; 10 |]
    in
    let y = AD.bmm x1 x2 in
    Algodiff.D.Maths.sum' y
  in
  let samples, directions = FD.generate_test_samples (1, 2000) n_samples in
  let threshold = 1E-4 in
  let eps = 1E-5 in
  let directions = Owl.Stats.shuffle directions in
  let directions = Array.sub directions 0 10 in
  FD.Reverse.check ~threshold ~order:`fourth ~eps ~directions ~f:ff samples |> fst
