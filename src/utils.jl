#

@noinline function _load_optdep(pkgid, msg)
    try
        Base.require(pkgid)
    catch err
        throw(ArgumentError("Error importing $pkgid for $msg: $err"))
    end
    return
end

@inline function load_optdep(pkgid, flag, msg)
    if !flag[]
        _load_optdep(pkgid, msg)
        flag[] = true
    end
    return
end
