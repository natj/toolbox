#system & io tools


# WARNING: nothing here works yet!!
# Everything you see is a work in progress

export ReadConf

#sreaddlm
#swritedlm

#print_shortest
#print_coolest

#ReadConf
############################################
#Read parameters from configure file
#
#input: file=path_to_configure_file
#list of variables to search eg. "teff"
#output: list of values read from conf-file
############################################
function ReadConf(file, vars...)
    f_stream=open(file)
    conf_values=readlines(f_stream)
    close(f_stream)
    var_arr=Any[]
    nvars=length(vars)

    for var in vars
        for line in conf_values
            #Reading all non-comment lines
            if !(ismatch(r"^\s*(?:#|$)",line))
                 m=search(line,var*"=")
                 if m != (0:-1);
                     var_s=strip(line[m[end]+1:end])
                     m2=search(var_s,".")
                     if m2 == (0:-1)
                         var_f=int64(var_s)
                     else
                         var_f=float64(var_s)
                     end
#                     println("found ",var," = ",var_f)
                     if nvars == 1
                         return var_f
                     else
                         push!(var_arr,var_f)
                     end
                 end
            end
        end
    end
    return var_arr
end


#function print_coolest(io::IO, x::FloatingPoint)
function print_coolest(io::IO, x::AbstractFloat)
    if isnan(x); return write(io, "NaN"); end
    if x < 0 write(io,'-') end
    if isinf(x); return write(io, "Inf"); end

    #smallest change
    #<100 use f
    #>100 use E
    #integer?

end
#print_coolest(io::IO, x::Union(FloatingPoint, Integer)) = print_coolest(io, float(x))

function swritedlm(io::IO, a::AbstractVecOrMat, dlm)
    pb = PipeBuffer()
    nr = size(a,1)
    nc = size(a,2)
    for i = 1:nr
        for j = 1:nc
            writedlm_cell(pb, a[i,j], dlm)
            j == nc ? write(pb,'\n') : print(pb,dlm)
        end
        (nb_available(pb) > (16*1024)) && write(io, takebuf_array(pb))
    end
    write(io, takebuf_array(pb))
    nothing
end
writedlm{T}(io::IO, a::AbstractArray{T,0}, dlm) = writedlm(io, reshape(a,1), dlm)


#function fwritedlm(io::IO, a::AbstractVecOrMat; header::ASCIIString="")
function fwritedlm(io::IO, a::AbstractVecOrMat; header::String="")
#    pb = PibeBuffer()
    nr = size(a,1)
    nc = size(a,2)

    hnts = split(header, " ")
    if lenght(types) != nc
        hnts = split(header, ",")
        if length(types) != nc
            error("Could not parse header, use space or comma as delimiter")
        end
    end

    headers = Array(String, nc)
    types = Array(String, nc)
    for i = 1:nc
        hnt = split(hnts[i], "%")
        headers[i] = hnt[1]
        types[i] = length(hnt) == 2 ? hnt[2] : ""
    end

    println(headers)
    println(types)

#    write(io, takebuf_array(pb))
#    nothing
end

