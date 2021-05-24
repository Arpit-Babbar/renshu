<TeXmacs|1.99.19>

<style|<tuple|tmarticle|old-lengths>>

<\body>
  <\hide-preamble>
    #Partial Derivative symbol

    <assign|pdv|<\macro|num|den>
      <frac|\<partial\> <arg|num>|\<partial\> <arg|den>>
    </macro>>

    <assign|pdv*|<\macro|degree|num|den>
      <frac|<math|\<partial\><rsup|<arg|degree>><arg|num>>|\<partial\>
      <math|<arg|den><rsup|<arg|degree>>>>
    </macro>>

    #Ordinary derivative symbol

    <assign|dv|<\macro|num|den>
      <frac|\<mathd\><arg|num>|\<mathd\><arg|den>>
    </macro>>

    <assign|dv*|<\macro|degree|num|den>
      <frac|<math|d<rsup|<arg|degree>><arg|num>>|d
      <math|<arg|den><rsup|<arg|degree>>>>
    </macro>>

    #Bold symbols - Start

    <assign|bC|<macro|\<b-C\>>>

    <assign|bd|<macro|\<b-d\>>>

    <assign|bD|<macro|\<b-D\>>>

    <assign|be|<macro|\<b-e\>>>

    <assign|bE|<macro|\<b-E\>>>

    <assign|bof|<macro|\<b-f\>>>

    <assign|bF|<macro|\<b-F\>>>

    <assign|bg|<macro|\<b-g\>>>

    <assign|bG|<macro|\<b-G\>>>

    <assign|bl|<macro|\<b-l\>>>

    <assign|bM|<macro|\<b-M\>>>

    <assign|bp|<macro|\<b-p\>>>

    <assign|bP|<macro|\<b-P\>>>

    <assign|bq|<macro|\<b-q\>>>

    <assign|bS|<macro|\<b-S\>>>

    <assign|br|<macro|\<b-r\>>>

    <assign|bu|<macro|\<b-u\>>>

    <assign|bv|<macro|\<b-v\>>>

    <assign|bw|<macro|\<b-w\>>>

    <assign|bx|<macro|\<b-x\>>>

    #Bold symbols-end
  </hide-preamble>

  <chapter*|Implementing Lax-Fredrich>

  Our eventual goal will be to solve the Riemann problem for Euler equations
  defined here

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <with|font-series|bold|Euler equations<em|<strong|<verbatim|>>>>

      Assume that the velocity is of the form <math|<around*|(|u,0,0|)>> and
      everything depends only on <math|x,t>. then the Euler equations in 1-D
      are given by

      <\equation*>
        U<rsub|t>+F<around*|(|U|)><rsub|x>=0
      </equation*>

      where

      <\equation*>
        U=<bmatrix|<tformat|<table|<row|<cell|\<rho\>>>|<row|<cell|\<rho\>u>>|<row|<cell|E>>>>>,<space|2em>F<around*|(|U|)>=<bmatrix|<tformat|<table|<row|<cell|\<rho\>u>>|<row|<cell|p+\<rho\>u<rsup|2>>>|<row|<cell|<around*|(|E+p|)>u>>>>>
      </equation*>

      and

      <\align>
        <tformat|<table|<row|<cell|\<rho\>=density,<space|1em>u=velocity,<space|1em>p=pressure>|<cell|>>|<row|<cell|E=total
        energy per unit volume=\<rho\>e+<frac|1|2>\<rho\>u<rsup|2>>|<cell|>>|<row|<cell|\<rho\>e=internal
        energy per unit volume>|<cell|>>|<row|<cell|e=internal energy per
        unit mass>|<cell|>>>>
      </align>

      The pressure <math|p> is related to the internal energy <math|e> by the
      caloric equation of state (????) <math|p=p<around*|(|\<rho\>,e|)>>; for
      a calorically ideal gas, <math|p=<around*|(|\<gamma\>-1|)>\<rho\>e>, so
      that

      <\equation*>
        p=<around*|(|\<gamma\>-1|)><around*|[|E-<frac|1|2>\<rho\>u<rsup|2>|]>.
      </equation*>
    </cell>>>>
  </wide-block>

  First, we will make a solver of Riemann problems for

  <\equation*>
    U<rsub|t>+A U<rsub|x>=0
  </equation*>

  where <math|U=<around*|(|U<rsub|1>,U<rsub|2>,U<rsub|3>|)>>,
  <math|U=U<around*|(|x,t|)>> and <math|A> is a diagonalizable matrix using
  Upwind and Lax-Friedrichs. First, let's do the Lax-Friedrichs flux.

  <subsection*|Lax-Friedrich>

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <subsubsection*|Dictionaries>

      <with|font-series|bold|<verbatim|equation>><space|1em>-<space|2em><verbatim|eq,
      flux, speed, compute_exact_soln!, name>

      <with|font-series|bold|<verbatim|problem>> - <verbatim|domain, nvar,
      initial_value, boundary_value, boundary_condition, final_time>

      <with|font-series|bold|<verbatim|parameter>><space|1em>-<space|2em><verbatim|grid_size,
      cfl, save_time_interval>

      <with|font-series|bold|<verbatim|scheme<with|font-series|bold|>>><space|2em>-<space|3em><verbatim|numerical_flux>

      \;

      <verbatim|nx/grid_size>
    </cell>>>>
  </wide-block>

  \;

  All functions in code will take arguments through dictionaries, while we
  will sometimes be explicit about the arguments here

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <subsubsection*|Grid file>

      <verbatim|make_grid(xmin, xmax, grid_size)> function takes these
      arguments to make a <verbatim|struct> of uniform grid of the interval
      [<verbatim|xmin,xmax>].
    </cell>>>>
  </wide-block>

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <subsubsection*|Run File for Linear Advection>

      <with|font-series|bold|User Inputs>

      <with|font-series|bold|<verbatim|equation>><space|1em>-
      <verbatim|fprime>

      <with|font-series|bold|<verbatim|problem>> - <verbatim|domain, nvar,
      initial_value, boundary_value, boundary_condition, final_time>

      <with|font-series|bold|<verbatim|parameter>><space|1em>-<space|2em><verbatim|grid_size,
      cfl, save_time_interval>

      <with|font-series|bold|<verbatim|scheme<with|font-series|bold|>>><space|2em>-<space|3em><verbatim|numerical_flux>

      \;

      The dictionary out of <with|font-series|bold|<verbatim|equation<with|font-series|bold|<verbatim|>>>>
      contains the <verbatim|struct> <verbatim|LinAdv>(classifying everything
      by classifying the PDE), physical <verbatim|flux>, advection
      <verbatim|speed> and a <verbatim|compute_exact_soln!> function.

      \;

      This is just the run file that is meant precisely for the PDE

      <\equation*>
        U<rsub|t>+A U<rsub|x>=0
      </equation*>

      So, here, we can specify what the matrix <math|A> will be and the file
      <verbatim|EqLinAdv.jl> will use the matrix to build the physical flux
      and numerical flux.
    </cell>>>>
  </wide-block>

  \;

  \;

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <\with|font-series|regular>
        <tabular*|<tformat|<table|<row|<cell|<with|font-series|bold|Main
        Algorithm> <verbatim|Solve>>>>>>
      </with>

      <\wide-block>
        <tformat|<table|<row|<\cell>
          <verbatim|compute lam, dt>

          For the general physical flux <math|F>, the <math|\<lambda\>> for
          Lax-Friedrich flux at a time level <verbatim|n> is chosen as

          <\equation*>
            \<lambda\>=max<rsub|j> \<sigma\><around*|(|F<rprime|'>*<around*|(|U<rsub|j><rsup|n>|)>|)>.
          </equation*>

          For this, we'd have to loop over all cells at each time step.

          And, then, we use it to compute the time step with safety factor
          CFL as

          <\equation*>
            d t=CFL<frac|\<mathLaplace\>x|\<lambda\>>,<space|1em>CFL\<leq\>1.
          </equation*>

          (Notice that <math|\<lambda\>> weirdly equals <math|<frac|d
          t|\<mathLaplace\>x>> in our case)

          \;
        </cell>>>>
      </wide-block>

      as
    </cell>>>>
  </wide-block>

  \;

  \;

  \;

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <subsubsection*|Exact Solution, Error evaluation for Linear System>

      If we are solving <math|U<rsub|t>+A U<rsub|x>=0>,
      <math|U<around*|(|x,0|)>=U<rsub|0><around*|(|x|)>><space|1em>for
      <math|\<Lambda\>=diag<around*|(|\<lambda\><rsub|1>,\<lambda\><rsub|2>,\<lambda\><rsub|3>|)>>
      with <math|\<lambda\><rsub|1>\<leq\>\<lambda\><rsub|2>\<leq\>\<lambda\><rsub|3>>
      such that we have <math|R> for which

      <\equation*>
        \<Lambda\>=R<rsup|-1>A R.
      </equation*>

      Then, if <math|U> is the exact solution to above equation,
      <math|W=R<rsup|-1>U> satisfies

      <\equation*>
        <tabular*|<tformat|<cwith|1|1|1|1|cell-halign|r>|<cwith|2|2|1|1|cell-halign|r>|<cwith|1|1|3|3|cell-halign|l>|<table|<row|<cell|<pdv|W<rsub|i>|t>+\<lambda\><rsub|i><pdv|W<rsub|i>|x>>|<cell|=>|<cell|0>|<cell|<space|1em>>|<cell|i=1,2,\<ldots\>,m>>|<row|<cell|W<rsub|i><around*|(|x,0|)>>|<cell|=>|<cell|<around*|(|R<rsup|-1>U<rsub|0><around*|(|x|)>|)><rsub|i>.>|<cell|>|<cell|>>>>>
      </equation*>

      So, we can easily compute <math|W> so that as
      <math|W<rsub|i><around*|(|x,t|)>=<around*|(|R<rsup|-1>U<rsub|0><around*|(|x-\<lambda\><rsub|i>t|)>|)><rsub|i>>

      Using this, we can plot the exact solution and compute convergence
      rate!\ 
    </cell>>>>
  </wide-block>

  \;

  \;

  \;

  <\wide-block>
    <tformat|<table|<row|<\cell>
      <with|font-series|bold|Improvements>

      Maybe, <verbatim|nvar> should be put in the equation file. Why will the
      user choose that? Of course, it doesn't matter for us because
      <verbatim|nvar> is just going to be 3 for this discussion.
    </cell>>>>
  </wide-block>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
    <associate|page-packet|1>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|?|1>>
    <associate|auto-2|<tuple|?|1>>
    <associate|auto-3|<tuple|?|1>>
    <associate|auto-4|<tuple|?|1>>
    <associate|auto-5|<tuple|?|2>>
    <associate|auto-6|<tuple|?|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|2fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-size|<quote|1.19>|Implementing
      Lax-Fredrich> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|1fn>

      <with|par-left|<quote|1tab>|Lax-Friedrich
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|2tab>|Grid file
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|Run File
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|Exact Solution, Error evaluation for Linear
      System <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|Variable Names
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>