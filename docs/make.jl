push!(LOAD_PATH,"../src/") 
using SeaGap
using Documenter

DocMeta.setdocmeta!(SeaGap, :DocTestSetup, :(using SeaGap); recursive=true)

makedocs(;
    modules=[SeaGap],
    authors="Fumiaki Tomita <39943988+f-tommy@users.noreply.github.com> and contributors",
    repo="https://github.com/f-tommy/SeaGap.jl/blob/{commit}{path}#{line}",
    sitename="SeaGap.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://f-tommy.github.io/SeaGap.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Dataformat"=> "Dataformat/intro.md",
        "Methodology" => [
            "Overview" => "Methodology/overview.md",
            "Travel-time calculation" => "Methodology/traveltime.md",
            "Kinematic array positioning" => "Methodology/kinematic.md",
            "Static array positioning" => "Methodology/static.md",
            "Individual transponder positioning" => "Methodology/single.md",
            "Static array positioning with sound speed gradients" => "Methodology/mcmcpvg.md",
          ],
        "Tutorials" => [
            "Forward calculation" => "Tutorials/forward.md",
            "Outlier removal" => "Tutorials/denoise.md",
            "Kinematic array positioning" => "Tutorials/tkinematic.md",
            "Static array positioning" => "Tutorials/tstatic.md",
            "Individual transponder positioning" => "Tutorials/tsingle.md",
            "MCMC array positioning" => "Tutorials/tmcmcpvg.md",
            "Time-series analysis" => "Tutorials/timeseries.md",
          ],
        "Others" => "Others/intro.md",
        "API" => "API.md"
    ],
)

deploydocs(;
    repo="github.com/f-tommy/SeaGap.jl",
    devbranch="main",
)
