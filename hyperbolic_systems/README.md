This is my first GPU code. It is a solver for 1-D Euler's equations written under the mentorship of [Julian Samaroo](https://github.com/jpsamaroo) during Hackathon-[JuliaCon 2022](https://juliacon.org/2022/).

I hope for it to serve as an introduction to GPU programming for CFD codes. Please raise an issue or start a discussion if you have any suggestions on how this repository can do so better.

I suggest doing the below exercise for hands-on practice. I wish to add some independent exercises at some point.

## Exercise 1

Put the variable `t` in a 1-element vector

## Exercise 2

The `compute_dt` function uses `@atomic` and thus does only one thread at a time. See the [JuliaCon 2021 tutorial](https://github.com/maleadt/juliacon21-gpu_workshop/blob/main/deep_dive/CUDA.ipynb) and figure a way to bypass this.

## Exercise 3

Implement second order FVM on GPU.
