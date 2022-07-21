variable2string(v) = String(Symbol(v))

function sort_symbols(symbol_array)
   string_array = collect(String.(symbol_array))
   sort!(string_array)
   return Symbol.(string_array)
end

function remove!(a, item)
   deleteat!(a, findall(x->x==item, a))
   nothing
end

function trim_fields!(fields, scheme::Scheme{<:Any, <:Any, <:Any, <:Any})
   remove!(fields, :limiter)
end

function trim_fields!(fields, param::Parameters{<:Any})
   nothing
end

function fieldnames_sorted(object)
   fields_unsorted = fieldnames(typeof(object))
   fields = sort_symbols(fields_unsorted)
   trim_fields!(fields, object)
   return fields
end

function get_sorted_field_values(object)
   fields = fieldnames_sorted(object)
   collect(  variable2string(getfield(object, field)) for field in fields  )
end

function get_filename(scheme, param)
   scheme_info = get_sorted_field_values(scheme)
   param_info = get_sorted_field_values(param)
   if "lwfr" in scheme_info
      remove!(x) = remove!(param_info, x)
      remove!.(["by degree", "RK11", "SSPRK22", "SSPRK33", "SSPRK54", "Tsit5",
                "RK4"])
   end
end
