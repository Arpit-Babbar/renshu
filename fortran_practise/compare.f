c This program compares two numbers, and introduces double precision #
        program compare
        implicit none 
c implicit none means that you have to specify the type
        real*8 a,b
c This is one way to initialize double precision numbers
        print*,'Enter the first number (a)'
        read*,a
        print*,'Enter second number (b)'
        read*,b
        if (a.lt.b)then
          print*,'a is less than b'
        else if(a.gt.b)then
          print*,'a is greater than b'
        else
          print*,'a is equal to b'
        endif
        end
      
