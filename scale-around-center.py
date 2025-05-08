import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
from vtkmodules.vtkFiltersCore import (
    vtkAppendFilter,
    vtkConnectivityFilter,
    vtkPolyDataConnectivityFilter,
    vtkDelaunay3D
)
import vtkmodules.vtkInteractionStyle
import vtk
import os
import sys

args = sys.argv[1:]

colors = vtkNamedColors()

filename = args[0]
reader = vtk.vtkSTLReader()
reader.SetFileName(filename)
reader.Update()

connectivityFilter = vtkPolyDataConnectivityFilter()
connectivityFilter.SetInputConnection(reader.GetOutputPort())
connectivityFilter.SetExtractionModeToSpecifiedRegions()
connectivityFilter.ColorRegionsOn()
components = list()

scale_factor = 5
idx = 0
while True:
    connectivityFilter.AddSpecifiedRegion(idx)
    connectivityFilter.Update()

    component = vtk.vtkPolyData()
    component.DeepCopy(connectivityFilter.GetOutput())

    data = connectivityFilter.GetOutput()
    print(f"{data=}")
    if component.GetNumberOfCells() <= 0:
        break
    
    components.append(component)
    connectivityFilter.DeleteSpecifiedRegion(idx)
    idx+=1

print(f"{len(components)=}")
if len(components) > 1:
    appendFilter = vtkAppendFilter()
    for i in range(len(components)):
        bounds = components[i].GetBounds()
        center = ((bounds[0] + bounds[1])/2,
                (bounds[2] + bounds[3])/2,
                (bounds[4] + bounds[5])/2)

        transform_to_center = vtkTransform()
        transform_to_center.Translate(-center[0], -center[1], -center[2])

        transformTranslateFilter = vtk.vtkTransformPolyDataFilter()
        transformTranslateFilter.SetInputConnection(reader.GetOutputPort())
        transformTranslateFilter.SetTransform(transform_to_center)
        transformTranslateFilter.Update()

        scale = vtkTransform()
        scale.Scale(scale_factor, scale_factor, scale_factor)

        transformScaleFilter = vtk.vtkTransformPolyDataFilter()
        transformScaleFilter.SetInputConnection(transformTranslateFilter.GetOutputPort())
        transformScaleFilter.SetTransform(scale)
        transformScaleFilter.Update()

        transform_from_center = vtkTransform()
        transform_from_center.Translate(center[0], center[1], center[2])

        transformScaleTransformedFilter = vtk.vtkTransformPolyDataFilter()
        transformScaleTransformedFilter.SetInputConnection(transformScaleFilter.GetOutputPort())
        transformScaleTransformedFilter.SetTransform(transform_from_center)
        transformScaleTransformedFilter.Update()

        appendFilter.AddInputData(transformScaleTransformedFilter.GetOutput())

    appendFilter.Update()

    transformed = appendFilter
    print(transformed)

else:
    for i in range(len(components)):
        bounds = components[i].GetBounds()
        center = ((bounds[0] + bounds[1])/2,
                (bounds[2] + bounds[3])/2,
                (bounds[4] + bounds[5])/2)

        transform_to_center = vtkTransform()
        transform_to_center.Translate(-center[0], -center[1], -center[2])

        transformTranslateFilter = vtk.vtkTransformPolyDataFilter()
        transformTranslateFilter.SetInputConnection(reader.GetOutputPort())
        transformTranslateFilter.SetTransform(transform_to_center)
        transformTranslateFilter.Update()

        scale = vtkTransform()
        scale.Scale(scale_factor, scale_factor, scale_factor)

        transformScaleFilter = vtk.vtkTransformPolyDataFilter()
        transformScaleFilter.SetInputConnection(transformTranslateFilter.GetOutputPort())
        transformScaleFilter.SetTransform(scale)
        transformScaleFilter.Update()

        transform_from_center = vtkTransform()
        transform_from_center.Translate(center[0], center[1], center[2])

        transformScaleTransformedFilter = vtk.vtkTransformPolyDataFilter()
        transformScaleTransformedFilter.SetInputConnection(transformScaleFilter.GetOutputPort())
        transformScaleTransformedFilter.SetTransform(transform_from_center)
        transformScaleTransformedFilter.Update()

    transformed = transformScaleTransformedFilter
    print(transformed)


if "-nopopup" not in args:
    originalMapper = vtkPolyDataMapper()
    originalMapper.SetInputConnection(reader.GetOutputPort())
    originalActor = vtkActor()
    originalActor.SetMapper(originalMapper)
    originalActor.GetProperty().SetColor(colors.GetColor3d('Blue'))

    transformedMapper = vtkDataSetMapper()
    transformedMapper.SetInputConnection(transformed.GetOutputPort())
    transformedActor = vtkActor()
    transformedActor.GetProperty().SetOpacity(0.5)
    transformedActor.SetMapper(transformedMapper)
    transformedActor.GetProperty().SetColor(colors.GetColor3d('Red'))

    renderer = vtkRenderer()
    renderer.AddActor(originalActor)
    renderer.AddActor(transformedActor)
    renderer.SetBackground(colors.GetColor3d('Green'))

    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderWindow.SetWindowName('TransformPolyData')
    renderWindow.Render()



    renderWindowInteractor.Start()
    

stlWriter = vtk.vtkSTLWriter()
stlWriter.SetFileName("fibrose_expandida.stl")
stlWriter.SetInputConnection(transformed.GetOutputPort())
stlWriter.Write()