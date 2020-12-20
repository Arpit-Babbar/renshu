c This program sums the first n numbers.
c Written in vim
      program array_sum
      implicit none
      integer i,nmax,n
c We initialize a constant variable with value 100
      parameter(nmax=100)
c Initialize array to store numbers to be added, variable to store sum
      real     val(nmax),sum
c Check how many numbers to add
      print*,'How many numbers to add?'
      read*,n

c Check that n is a valid input
      if (n.gt.nmax)then
        print*,'Array size nmax is not sufficient.'
        print*,'Increase nmax and then run'
        stop
      endif

      if (n.lt.0)then
        print*,'Error: give a positive integer'
        stop
      endif
c Initialize the numbers to add
c Writing a do loop without increment size assumes increment=1
      do i=1,n
        val(i)=i
      enddo

c Add the numbers
        sum = 0.0
        do i=1,n
          sum = sum + val(i)
        enddo

        print*,'sum =',sum

        end program array_sum

