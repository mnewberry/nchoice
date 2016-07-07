
(* Use GSL's random number generator *)
let rng = Gsl.Rng.make (Gsl.Rng.default ())

let init seed = Gsl.Rng.set rng (Nativeint.of_int seed)
let int i = if i = 1 then 0 else Gsl.Rng.uniform_int rng i
let float scale = scale *. Gsl.Rng.uniform rng
let lcalphastr n = String.map (fun _ -> Char.chr (int 26 + 97)) 
  (String.make n ' ')

let binom n p = Gsl.Randist.binomial rng p n
let multinom n ps = Array.to_list 
  (Gsl.Randist.multinomial rng n (Array.of_list ps))

(* Use ocaml's random number generator
let init = Random.init
let int i = if i = 1 then 0 else Random.int i
let float = Random.float  *)

let exp lambda = log (float 1.0) /. (0. -. lambda) 
