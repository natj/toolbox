# Small tools for debugging

export catch_NaN


#nan checker
function catch_NaN(M)
    l = length(M)
    for i = 1:l
        if isnan(M[i])
            error("Found NaN in $i")
        end
        if M[i] < 0.0
            error("Found negative value at $i")
        end
#        if M[i] == 0.0
#            warn("Found zero value at $i")
#        end
    end
    return false
end