open Owl
open Bmo

let n = 5

let test () =
  let x =
    let x = Arr.gaussian [| 5; 3; n; n |] in
    D.print_dim x;
    let xt = Arr.(transpose ~axis:[| 0; 1; 3; 2 |] x) in
    D.print_dim xt;
    let x = D.bmm x xt in
    let e = Mat.eye n in
    Arr.(x + e)
  in
  let z = D.bchol ~upper:false x in
  ignore z


let test2 () =
  let e = Mat.eye n in
  let x =
    Array.init 100 (fun _ ->
        let z = Mat.(gaussian n n) in
        Mat.(e + (z *@ transpose z)))
  in
  let x' = Arr.stack ~axis:0 x in
  let check ~upper =
    let z = Array.map (Linalg.D.chol ~upper) x |> Arr.stack ~axis:0 in
    let z' = D.bchol ~upper x' in
    assert (Arr.(approx_equal z z'))
  in
  check ~upper:true;
  check ~upper:false


let () =
  test ();
  test2 ()
