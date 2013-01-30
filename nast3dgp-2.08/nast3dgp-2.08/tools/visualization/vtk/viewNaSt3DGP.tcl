catch {load vtktcl}

set DIRECTORY $argv

set FLAG_FILE $DIRECTORY.vtk.fl
set PRESSURE_FILE $DIRECTORY.vtk.p

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
    ren1 SetBackground 1 1 1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin
vtkRectilinearGridReader flagReader
    flagReader SetFileName $FLAG_FILE


if { [file exists $FLAG_FILE] == 1 } { 
    flagReader SetFileName $FLAG_FILE
    flagReader Update
} else {
    return
}
set dims [[flagReader GetOutput] GetDimensions]
set dimx [lindex $dims 0]
set dimy [lindex $dims 1]
set dimz [lindex $dims 2]

# renWin SetSize [expr $dimx * 8] [expr $dimy * 8]

vtkOutlineFilter outlineFilter
vtkOutlineSource outlineSource
vtkPolyDataMapper outlineMapper
vtkActor outline 

outlineFilter SetInput [flagReader GetOutput]
outlineMapper SetInput [outlineFilter GetOutput]

[outline GetProperty] SetAmbient 1.0
[outline GetProperty] SetDiffuse 1         
[outline GetProperty] SetColor 0.1 0.1 0.1

outline SetMapper outlineMapper
ren1 AddActor outline

#
# IsoFlaeche des Flagfeldes 
#
vtkContourFilter flagIso   
vtkPolyDataMapper flagIsoMapper
vtkActor flagIsoActor
  flagIso SetInput [flagReader GetOutput]
  flagIso SetValue 0 1  

vtkPolyDataNormals normals
    normals SetInput [flagIso GetOutput]
    normals SetFeatureAngle  45
    flagIsoMapper SetInput [normals GetOutput]
    flagIsoMapper ScalarVisibilityOff
  
    [flagIsoActor GetProperty] SetDiffuse 1.0
    [flagIsoActor GetProperty] SetColor 0.9 0.9 0.9
    [flagIsoActor GetProperty] SetOpacity 1
    flagIsoActor SetMapper flagIsoMapper
 
    ren1 AddActor flagIsoActor


#*********************************************
#
#       Schnitt-Ebenen des Feldes 
#
#*********************************************


vtkRectilinearGridGeometryFilter slice1        
vtkRectilinearGridGeometryFilter slice2       
vtkRectilinearGridGeometryFilter slice3        
vtkRectilinearGridGeometryFilter slice4   
vtkPolyDataMapper sliceMapper1  
vtkPolyDataMapper sliceMapper2 
vtkPolyDataMapper sliceMapper3     
vtkPolyDataMapper sliceMapper4   
vtkActor sliceActor1  
vtkActor sliceActor2
vtkActor sliceActor3   
vtkActor sliceActor4      
vtkContourFilter cont2
vtkContourFilter cont3
vtkContourFilter cont4

  # Lesen des u-Feldes aus dem .u File ->slice
  
vtkRectilinearGridReader u1      
vtkRectilinearGridReader u2       
vtkLookupTable lut    
vtkVRMLExporter exp  

    u1 SetFileName $FLAG_FILE          
    u2 SetFileName $PRESSURE_FILE

    slice1 SetInput [u1 GetOutput]        
    slice2 SetInput [u2 GetOutput]     
    slice3 SetInput [u2 GetOutput]     
    slice4 SetInput [u2 GetOutput]   
       
# ISO-Flaeche fuer U, V, W, P, T oder C
vtkContourFilter iso1   
vtkContourFilter iso2    
    iso1 SetInput [u1 GetOutput] 
    iso1 SetValue 0 -.2         
    iso2 SetInput [u1 GetOutput] 
    iso2 SetValue 0 0.0     

vtkPolyDataNormals normals1    
vtkPolyDataNormals normals2     
    normals1 SetInput [iso1 GetOutput]
    normals1 SetFeatureAngle  45
    normals2 SetInput [iso2 GetOutput]
    normals2 SetFeatureAngle  45 

vtkPolyDataMapper isoMapper1
vtkPolyDataMapper isoMapper2 
    isoMapper1 SetInput [normals1 GetOutput]
    isoMapper1 SetLookupTable lut
    isoMapper1 ScalarVisibilityOff
    isoMapper1 SetScalarRange -0.5 0.5
    isoMapper2 SetInput [normals2 GetOutput]
    isoMapper2 SetLookupTable lut 
    isoMapper2 ScalarVisibilityOff
    isoMapper2 SetScalarRange -0.5 0.5 
 
    lut SetHueRange 0.667 0.0

vtkActor isoActor1   
vtkActor isoActor2  

    isoActor1 SetMapper isoMapper1
    isoActor2 SetMapper isoMapper2

    [isoActor1 GetProperty] SetColor 0 0 1
    [isoActor1 GetProperty] SetAmbient -0.2
    [isoActor1 GetProperty] SetOpacity 1
    [isoActor2 GetProperty] SetColor 0 0 1  
    [isoActor2 GetProperty] SetAmbient 0.0
    [isoActor2 GetProperty] SetOpacity 1  
 
    slice2 SetExtent 0 0 0 0 0 0
    slice3 SetExtent 0 0 0 0 0 0  
    slice4 SetExtent 0 0 0 0 0 0  

    sliceMapper2 SetLookupTable lut 
    sliceMapper3 SetLookupTable lut   
    sliceMapper4 SetLookupTable lut

    sliceMapper2 SetInput [slice2 GetOutput]
    sliceMapper3 SetInput [slice3 GetOutput]
    sliceMapper4 SetInput [slice4 GetOutput]  
    sliceMapper2 ScalarVisibilityOn
    sliceMapper3 ScalarVisibilityOn 
    sliceMapper4 ScalarVisibilityOn 

    sliceMapper2 SetScalarRange -.2 .2
    sliceMapper3 SetScalarRange -.2 .2
    sliceMapper4 SetScalarRange -.2 .2 

    lut SetHueRange 0.667 0.0

    sliceActor2 SetMapper sliceMapper2
    sliceActor3 SetMapper sliceMapper3  
    sliceActor4 SetMapper sliceMapper4

    ren1 AddActor sliceActor2
    ren1 AddActor sliceActor3   
    ren1 AddActor sliceActor4 

    ren1 AddActor isoActor2 


# render the image
#
# iren SetUserMethod {wm deiconify .vtkInteract}
set cam1 [ren1 GetActiveCamera]
# $cam1 Zoom 1.5
iren Initialize

#
# Create user interface
#
frame .mbar -borderwidth 1 -relief raised
pack .mbar -fill x

vtkPNMWriter PNMwriter
vtkWindowToImageFilter win2iF
set nr [expr 1]
button .mbar.m -text "Exit" -state normal -command exit
button .mbar.photo -text "ppm" -state normal -command {
    u1 Update
    renWin Render
    win2iF SetInput renWin
    win2iF Modified
    PNMwriter SetInput [win2iF GetOutput]
    PNMwriter SetFileName snapshot_$nr.ppm
    PNMwriter Write
    incr nr
}


menubutton .mbar.freesurf -text "IsoSurf" -menu .mbar.freesurf.m
menu .mbar.freesurf.m
.mbar.freesurf.m add command -label "U" -command {  
    if { [file exists $DIRECTORY.vtk.u] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.u
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 1 0 0
	renWin Render
    } else {
        return
    }
}
.mbar.freesurf.m add command -label "V" -command {  
    if { [file exists $DIRECTORY.vtk.v] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.v
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 0 1 0
        renWin Render
    } else {
        return
    }
}
.mbar.freesurf.m add command -label "W" -command {  
    if { [file exists $DIRECTORY.vtk.w] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.w
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 1 1 0
        renWin Render
    } else {
        return
    }
}
.mbar.freesurf.m add command -label "P" -command {  
    if { [file exists $DIRECTORY.vtk.p] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.p
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 1 0 1
        renWin Render
    } else {
        return
    }
}
.mbar.freesurf.m add command -label "Tp" -command {  
    if { [file exists $DIRECTORY.vtk.t] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.t
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 0 1 1
        renWin Render
    } else {
        return
    }
}
.mbar.freesurf.m add command -label "Ch" -command {  
    if { [file exists $DIRECTORY.vtk.c0] == 1 } {
	set op [[isoActor2 GetProperty] GetOpacity]
        u1 SetFileName $DIRECTORY.vtk.c0
	[isoActor2 GetProperty] SetOpacity $op
	[isoActor2 GetProperty] SetColor 1 0 0
        renWin Render
    } else {
        return
    }
}

pack .mbar.m .mbar.photo -side left
pack .mbar.freesurf -side right

menubutton .mbar.slice -text "Slice" -menu .mbar.slice.m
menu .mbar.slice.m
.mbar.slice.m add command -label "U" -command {
    if { [file exists $DIRECTORY.vtk.u] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.u
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
	slice3 SetExtent 0 $dimx $y $y 0 $dimz
	slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}
.mbar.slice.m add command -label "V" -command {
    if { [file exists $DIRECTORY.vtk.v] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.v
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
        slice3 SetExtent 0 $dimx $y $y 0 $dimz
        slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}
.mbar.slice.m add command -label "W" -command {
    if { [file exists $DIRECTORY.vtk.w] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.w
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
        slice3 SetExtent 0 $dimx $y $y 0 $dimz
        slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}
.mbar.slice.m add command -label "P" -command {
    if { [file exists $DIRECTORY.vtk.p] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.p
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
        slice3 SetExtent 0 $dimx $y $y 0 $dimz
        slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}
.mbar.slice.m add command -label "Tp" -command {
    if { [file exists $DIRECTORY.vtk.t] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.t
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
        slice3 SetExtent 0 $dimx $y $y 0 $dimz
        slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}
.mbar.slice.m add command -label "Ch" -command {
    if { [file exists $DIRECTORY.vtk.c0] == 1 } {
        u2 SetFileName $DIRECTORY.vtk.c0
	slice2 SetExtent 0 $dimx 0 $dimy $z $z
        slice3 SetExtent 0 $dimx $y $y 0 $dimz
        slice4 SetExtent $x $x 0 $dimy 0 $dimz
        renWin Render
    } else {
        return
    }
}

pack .mbar.m .mbar.photo -side left
pack .mbar.slice -side right
pack .mbar.freesurf -side right


set fra [expr 0]
set inv [expr 0]
frame .f2
label .f2.l2 -text "Surf Op."
scale .f2.op -from 0 -to 100  \
    -orient horizontal -command SetOpacity
button .f2.1 -text "Frame" -width 3 -state normal -command {
    [outline GetProperty] SetOpacity $fra
    set fra [expr 1-$fra]
    renWin Render
}

frame .f3
label .f3.l3 -text "Flag Op."
scale .f3.opflag -from 0 -to 100  \
    -orient horizontal -command SetOpacityFlag
button .f3.1 -text " Inv " -width 3 -state normal -command {
    [outline GetProperty] SetColor [expr 1-$inv] [expr 1-$inv] [expr 1-$inv]
    ren1 SetBackground $inv $inv $inv
    set inv [expr 1-$inv]
    renWin Render
}

frame .f4
label .f4.tolab -text "IsoSurfValue:" 
entry .f4.to -width 10
button .f4.1 -text "On/Off" -width 3 -state normal -command {
    set op [[isoActor2 GetProperty] GetOpacity]
    if { $op == 0 } {set op [expr [.f2.op get]/100.0]} else {set op [expr 0]}
    [isoActor2 GetProperty] SetOpacity $op
    renWin Render
}

frame .fr
label .fr.mintolab -text " Range min:" 
entry .fr.minto -width 18
 
frame .fl 
label .fl.maxtolab -text "Range max:" 
entry .fl.maxto -width 18

set z [expr 0]
set onz [expr 1]
frame .f3b
button .f3b.1 -text "SliceZ +" -state normal -command {
    if { $z >= $dimz-1 } { } else { incr z }
    slice2 SetExtent 0 $dimx 0 $dimy $z $z
    [sliceActor2 GetProperty] SetOpacity $onz
    renWin Render
}
button .f3b.2 -text "SliceZ -" -state normal -command {
    if { $z <= 0 } { } else { set z [expr $z-1] } 
    slice2 SetExtent 0 $dimx 0 $dimy $z $z
    [sliceActor2 GetProperty] SetOpacity $onz
    renWin Render
}
button .f3b.3 -text "On/Off" -width 4 -state normal -command {
    sliceMapper2 SetInput [slice2 GetOutput]
    slice2 SetExtent 0 $dimx 0 $dimy $z $z
    if { $onz == 1 } {set onz [expr 0]} else {set onz [expr 1]}
    [sliceActor2 GetProperty] SetOpacity $onz
    renWin Render
}


set y [expr 0]
set ony [expr 1]
frame .f3c
button .f3c.1 -text "SliceY +" -state normal -command {
    if { $y >= $dimy-1 } { } else { incr y }
    slice3 SetExtent 0 $dimx $y $y 0 $dimz
    [sliceActor3 GetProperty] SetOpacity $ony
    renWin Render
}
button .f3c.2 -text "SliceY -" -state normal -command {
    if { $y <= 0 } { } else { set y [expr $y-1] }
    slice3 SetExtent 0 $dimx $y $y 0 $dimz
    [sliceActor3 GetProperty] SetOpacity $ony
    renWin Render
}
button .f3c.3 -text "On/Off" -width 4 -state normal -command {
    sliceMapper3 SetInput [slice3 GetOutput]
    slice3 SetExtent 0 $dimx $y $y 0 $dimz
    if { $ony == 1 } {set ony [expr 0]} else {set ony [expr 1]}
    [sliceActor3 GetProperty] SetOpacity $ony
    renWin Render
}


set x [expr 0]
set onx [expr 1]
frame .f3d
button .f3d.1 -text "SliceX +" -state normal -command {
    if { $x >= $dimx-1 } { } else { incr x }
    slice4 SetExtent $x $x  0 $dimy 0 $dimz
    [sliceActor4 GetProperty] SetOpacity $onx
    renWin Render
}
button .f3d.2 -text "SliceX -" -state normal -command {
    if { $x <= 0 } { } else { set x [expr $x-1] }
    slice4 SetExtent $x $x  0 $dimy 0 $dimz
    [sliceActor4 GetProperty] SetOpacity $onx
    renWin Render
}
button .f3d.3 -text "On/Off" -width 4 -state normal -command {
    sliceMapper4 SetInput [slice4 GetOutput]
    slice4 SetExtent $x $x  0 $dimy 0 $dimz
    if { $onx == 1 } {set onx [expr 0]} else {set onx [expr 1]}
    [sliceActor4 GetProperty] SetOpacity $onx
    renWin Render
}

set opacity [[isoActor2 GetProperty] GetOpacity]
.f2.op set [expr [lindex $opacity 0] * 100.0]

set opacityflag [[flagIsoActor GetProperty] GetOpacity]
.f3.opflag set [expr [lindex $opacityflag 0] * 100.0]

pack .f2.l2 .f2.op .f2.1 -side left
pack .f3.l3 .f3.opflag .f3.1 -side left
pack .f4.1 .f4.tolab .f4.to -side left
pack .fr.mintolab .fr.minto -side left
pack .fl.maxtolab .fl.maxto -side left
pack .f2 .f3 .f4 .f3d .f3c .f3b .fr .fl 
pack .f3b.1 .f3b.2 .f3b.3 -side left
pack .f3c.1 .f3c.2 .f3c.3 -side left
pack .f3d.1 .f3d.2 .f3d.3 -side left

bind .f4.to <KeyPress-Return> {
    iso2 SetValue 0 [.f4.to get]
    renWin Render
}

bind .fr.minto <KeyPress-Return> {
    sliceMapper2 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    sliceMapper3 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    sliceMapper4 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    renWin Render
}

bind .fl.maxto <KeyPress-Return> {
    sliceMapper2 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    sliceMapper3 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    sliceMapper4 SetScalarRange [expr [.fr.minto get]] [expr [.fl.maxto get]]
    renWin Render
}

proc SetOpacity {value} {
    [isoActor2 GetProperty] SetOpacity [expr [.f2.op get]/100.0]
    renWin Render
}

proc SetOpacityFlag {value} {
    [flagIsoActor GetProperty] SetOpacity [expr [.f3.opflag get]/100.0]
    renWin Render
}

