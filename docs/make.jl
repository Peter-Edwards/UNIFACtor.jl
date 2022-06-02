using Documenter
using UNIFACtor

makedocs(
    sitename = "UNIFACtor",
    format = Documenter.HTML(),
    modules = [UNIFACtor]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
