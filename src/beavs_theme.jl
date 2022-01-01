set_theme!(
    Theme(
        palette = (color=[c for c in ColorSchemes.Dark2_5], marker=[:circle, :utriangle, :cross, :rect, :diamond, :dtriangle, :pentagon, :xcross]),
        textcolor = :gray40,
        linewidth=4,
        fontsize=20,
        font="Courier",
        resolution = (520, 400),
        Axis = (
            backgroundcolor = RGB(0.96, 1.0, 0.98),
            xgridcolor = (:black, 0.15),
            ygridcolor = (:black, 0.15),
            leftspinevisible = false,
            rightspinevisible = false,
            ygridstyle=:dash,
            xgridstyle=:dash,
            bottomspinevisible = false,
            topspinevisible = false,
            xminorticksvisible = false,
            yminorticksvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            xlabelpadding = 3,
            ylabelpadding = 3
        ),
        Legend = (
            framevisible = true,
            titlehalign=:left,
            titlesize=16,
            labelsize=16,
            framecolor=(:black, 0.5)
            # padding = (1, 0, 0, 0),
        ),
        Axis3 = (
            xgridcolor = (:black, 0.07),
            ygridcolor = (:black, 0.07),
            zgridcolor = (:black, 0.07),
            xspinesvisible = false,
            yspinesvisible = false,
            zspinesvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
        ),
        Colorbar = (
            ticksvisible = false,
            spinewidth = 0,
            ticklabelpad = 5,
        )
    )
)
