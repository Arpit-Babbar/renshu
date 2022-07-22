using TimerOutputs
using LoopVectorization
using UnPack

function create_container(n)
   arr1, arr2, arr3 = [rand(n,n) for _=1:3]
   container = (; arr1, arr2, arr3 )
   return container
end

function add_stuff(n)
   container = create_container(n)
   @timeit_debug "This does nothing" begin
      @unpack arr1, arr2, arr3 = container
   end # timer
   # arr1, arr2, arr3 = [rand(n,n) for _=1:3]
   @turbo @. arr3 = arr1 + arr2
   return nothing
end

add_stuff(10)
