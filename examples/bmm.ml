open Owl
open Bmo

let test1 () =
  let nbatch = 100 in
  let m = 3 in
  let k = 2 in
  let n = 5 in
  let x = Array.init nbatch (fun _ -> Mat.gaussian m k) in
  let y = Array.init nbatch (fun _ -> Mat.gaussian k n) in
  let z = D.bmm Arr.(stack ~axis:0 x) Arr.(stack ~axis:0 y) in
  let z' = Array.(map2 (fun x y -> Mat.(x *@ y)) x y) |> Arr.stack ~axis:0 in
  Arr.print z;
  D.print_dim z;
  Arr.print z';
  D.print_dim z';
  assert (Arr.approx_equal ~eps:1E-6 z z')


let test2 () =
  let m = 3 in
  let n = 5 in
  let k = 2 in
  let x = Arr.gaussian [| 2; 3; 6; m; k |] in
  let y = Arr.gaussian [| 2; 3; 6; k; n |] in
  let z = D.bmm x y in
  D.print_dim z;
  let shp = Arr.shape z in
  assert (shp = [| 2; 3; 6; m; n |])


let () =
  test1 ();
  test2 ()
