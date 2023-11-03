function StrongWolfe1(fk, fj, gk, pk, alphaj; c1=1e-4)
    if fk - fj ≥ c1*alphaj*gk'*pk
        return true
    else
        return false
    end

    @error "floc"
end 

function StrongWolfe2(gk, gj, pk, alphaj; c2=0.9)
    if abs(gj'*pk) ≤ c2*abs(gk'*pk)
        return true
    else
        return false
    end

    @error "floc"
end