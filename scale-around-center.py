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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from pathlib import Path

# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
    VTK_DOUBLE_MAX,
    vtkPoints
)
from vtkmodules.vtkCommonCore import (
    VTK_VERSION_NUMBER,
    vtkVersion
)
from vtkmodules.vtkCommonDataModel import (
    vtkIterativeClosestPointTransform,
    vtkPolyData
)
from vtkmodules.vtkCommonTransforms import (
    vtkLandmarkTransform,
    vtkTransform
)
from vtkmodules.vtkFiltersGeneral import (
    vtkOBBTree,
    vtkTransformPolyDataFilter
)
from vtkmodules.vtkFiltersModeling import vtkHausdorffDistancePointSetFilter
from vtkmodules.vtkIOGeometry import (
    vtkBYUReader,
    vtkOBJReader,
    vtkSTLReader
)
from vtkmodules.vtkIOLegacy import (
    vtkPolyDataReader,
    vtkPolyDataWriter
    )
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.vtkInteractionWidgets import (
    vtkCameraOrientationWidget,
    vtkOrientationMarkerWidget
)
from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)

from vtkmodules.vtkFiltersHybrid import (
    vtkPCAAnalysisFilter,
    vtkProcrustesAlignmentFilter,
)
from vtkmodules.vtkFiltersGeneral import (
    vtkMultiBlockDataGroupFilter,
    vtkTransformPolyDataFilter,
)

from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

# def procrustes_aligment(source, target):
#     procrustes = vtkProcrustesAlignmentFilter()
#     group = vtkMultiBlockDataGroupFilter()
#     group.AddInputConnection(source)
#     group.AddInputConnection(target)
#     procrustes.SetInputConnection(group.GetOutputPort())
#     procrustes.GetLandmarkTransform().SetModeToRigidBody()
#     procrustes.Update()
#     return procrustes
def procrustes_aligment(source_polydata, target_polydata):
    group = vtkMultiBlockDataSet()
    group.SetNumberOfBlocks(2)
    group.SetBlock(0, source_polydata)
    group.SetBlock(1, target_polydata)

    procrustes = vtkProcrustesAlignmentFilter()
    procrustes.SetInputData(group)
    procrustes.GetLandmarkTransform().SetModeToRigidBody()
    procrustes.Update()

    return procrustes.GetOutput().GetBlock(0)


# def align_bounding_boxes(source, target):
#     # Use OBBTree to create an oriented bounding box for target and source
#     source_obb_tree = vtkOBBTree()
#     source_obb_tree.SetDataSet(source)
#     source_obb_tree.SetMaxLevel(1)
#     source_obb_tree.BuildLocator()

#     target_obb_tree = vtkOBBTree()
#     target_obb_tree.SetDataSet(target)
#     target_obb_tree.SetMaxLevel(1)
#     target_obb_tree.BuildLocator()

#     source_landmarks = vtkPolyData()
#     source_obb_tree.GenerateRepresentation(0, source_landmarks)

#     target_landmarks = vtkPolyData()
#     target_obb_tree.GenerateRepresentation(0, target_landmarks)

#     lm_transform = vtkLandmarkTransform()
#     lm_transform.SetModeToRigidBody()
#     # lm_transform.SetModeToSimilarity()
#     lm_transform.SetTargetLandmarks(target_landmarks.GetPoints())
#     best_distance = VTK_DOUBLE_MAX
#     best_points = vtkPoints()
#     best_distance = best_bounding_box(
#         "X",
#         target,
#         source,
#         target_landmarks,
#         source_landmarks,
#         best_distance,
#         best_points)
#     best_distance = best_bounding_box(
#         "Y",
#         target,
#         source,
#         target_landmarks,
#         source_landmarks,
#         best_distance,
#         best_points)
#     best_distance = best_bounding_box(
#         "Z",
#         target,
#         source,
#         target_landmarks,
#         source_landmarks,
#         best_distance,
#         best_points)

#     lm_transform.SetSourceLandmarks(best_points)
#     lm_transform.Modified()

#     lm_transform_pd = vtkTransformPolyDataFilter()
#     lm_transform_pd.SetInputData(source)
#     lm_transform_pd.SetTransform(lm_transform)
#     lm_transform_pd.Update()

#     source.DeepCopy(lm_transform_pd.GetOutput())

#     return


# def best_bounding_box(axis, target, source, target_landmarks, source_landmarks, best_distance, best_points):
#     distance = vtkHausdorffDistancePointSetFilter()
#     test_transform = vtkTransform()
#     test_transform_pd = vtkTransformPolyDataFilter()
#     lm_transform = vtkLandmarkTransform()
#     lm_transform_pd = vtkTransformPolyDataFilter()

#     # lm_transform.SetModeToRigidBody()
#     lm_transform.SetModeToSimilarity()
#     lm_transform.SetTargetLandmarks(target_landmarks.GetPoints())

#     source_center = source_landmarks.GetCenter()

#     delta = 90.0
#     for i in range(0, 4):
#         angle = delta * i
#         print(f"angle = {angle}")

#         # Rotate about center
#         test_transform.Identity()
#         test_transform.Translate(source_center[0], source_center[1], source_center[2])
#         if axis == "X":
#             test_transform.RotateX(angle)
#         elif axis == "Y":
#             test_transform.RotateY(angle)
#         else:
#             test_transform.RotateZ(angle)
#         test_transform.Translate(-source_center[0], -source_center[1], -source_center[2])

#         test_transform_pd.SetTransform(test_transform)
#         test_transform_pd.SetInputData(source_landmarks)
#         test_transform_pd.Update()

#         lm_transform.SetSourceLandmarks(test_transform_pd.GetOutput().GetPoints())
#         lm_transform.Modified()

#         lm_transform_pd.SetInputData(source)
#         lm_transform_pd.SetTransform(lm_transform)
#         lm_transform_pd.Update()

#         distance.SetInputData(0, target)
#         distance.SetInputData(1, lm_transform_pd.GetOutput())
#         distance.Update()

#         test_distance = distance.GetOutput(0).GetFieldData().GetArray("HausdorffDistance").GetComponent(0, 0)
#         if test_distance < best_distance:
#             best_distance = test_distance
#             best_points.DeepCopy(test_transform_pd.GetOutput().GetPoints())

#     return best_distance

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

scale_factor = 2
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
        transformTranslateFilter.SetInputData(components[i])
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
        procrustes_aligment(transformScaleTransformedFilter.GetOutput(), components[i])
        # align_bounding_boxes(transformScaleTransformedFilter.GetOutput(), components[i])

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

    scaled_polydata = transformScaleTransformedFilter.GetOutput()
    original_polydata = reader.GetOutput()

    aligned_polydata = procrustes_aligment(scaled_polydata, original_polydata)
    transformed = vtk.vtkPolyData()
    transformed.DeepCopy(aligned_polydata)
    # align_bounding_boxes(transformed.GetOutput(), reader.GetOutput())
    print(transformed)


if "-nopopup" not in args:

    print("Iniciando visualização VTK...")
    originalMapper = vtkPolyDataMapper()
    originalMapper.SetInputConnection(reader.GetOutputPort())
    originalActor = vtkActor()
    originalActor.SetMapper(originalMapper)
    originalActor.GetProperty().SetColor(colors.GetColor3d('Blue'))

    transformedMapper = vtkDataSetMapper()
    transformedMapper.SetInputData(transformed)
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