import Pkg


# tries to load a Python dependency; on failure, adds the dependency via Conda
function check_add_dep(pkg; channel="")
    try
        @info "Checking dependency: $pkg"
        pyimport(pkg)
    catch
        @info "Installing $pkg..."
        Conda.add(pkg, channel=channel)
        pyimport(pkg)
        @info "$pkg install verified."
    end
end


# check for a working Python environment; if none, install one
try
    @info "Checking PyCall."
    using PyCall
    pyimport("sys")
    @info "Python environment verified."
catch
    @info "Setting up Python environment..."
    ENV["PYTHON"]=""
    Pkg.add("PyCall")
    Pkg.build("PyCall")
    using PyCall
    pyimport("sys")
    @info "PyCall verified."
end

# check for Conda; if not found, install it
try
    @info "Checking Conda."
    using Conda
    @info "Conda verified."
catch
    @info "Installing Conda..."
    Pkg.add("Conda")
    using Conda
    @info "Conda verified."
end

# check for Controlz; if not found, add it
try
    @info "Checking Controlz."
    using Controlz
    @info "Controlz verified."
catch # if not, install it
    @info "Installing Controlz..."
    Pkg.add("Controlz")
    Pkg.build("Controlz")
    using Controlz
    @info "Controlz install verified."
end

check_add_dep("matplotlib")

@info "Setup complete!"
