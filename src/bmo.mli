open Owl

module D : sig
  type arr = Arr.arr

  val print_dim : arr -> unit
  val bmm : arr -> arr -> arr
  val bchol : ?upper:bool -> arr -> arr
end

module AD : sig
  type t = Owl.Algodiff.D.t

  val print_dim : t -> unit
  val bmm : t -> t -> t
end
