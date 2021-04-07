module 

using WriteVTK

#-------------------------------------------------------------------------------
# Return string of the form base_name00c with total number of digits = ndigits
#-------------------------------------------------------------------------------
function get_filename(base_name, ndigits, c)
   if c > 10^ndigits - 1
       println("get_filename: Not enough digits !!!")
       println("   ndigits =", ndigits)
       println("   c       =", c)
       exit()
   end
   number = lpad(c, ndigits, "0")
   return string(base_name, number)
end

#-------------------------------------------------------------------------------
function write_vtk(base_name, fcount, iter, time, grid, z, ndigits=3)
  # Clear and re-create output directory
  if fcount == 0
     run(`rm -rf output`)
     run(`mkdir output`)
  end
  filename = get_filename(base_name, ndigits, fcount)
  filename = string("output/", filename)
  vtk = vtk_grid(filename, grid.xc, grid.yc)
  nx, ny = grid.size
  vtk["Solution"] = @view z[1:nx,1:ny]
  vtk["CYCLE"] = iter
  vtk["TIME"] = time
  out = vtk_save(vtk)
  println("Wrote file ", out[1])
  fcount += 1
  return fcount
end
