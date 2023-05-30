# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonizerDockWidget
                                A QGIS plugin
 A plugin to help make Vision Zero Polygons
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                            -------------------
        begin                : 2023-05-27
        git sha              : $Format:%H$
        copyright            : (C) 2023 by DTS
        email                : frank.hereford@austintexas.gov
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
import json
#import sys
#sys.path.insert(0, '/Users/frank/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/polygonizer/')

from qgis.PyQt import QtGui, QtWidgets, uic
from qgis.PyQt.QtCore import pyqtSignal, QVariant

#from qgis.core import (
    #QgsProject,
    #QgsVectorLayer
#)

from qgis.core import *
from qgis.utils import *

#from qgis.utils import (
    #iface
#)

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'polygonizer_dockwidget_base.ui'))


class PolygonizerDockWidget(QtWidgets.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, parent=None):
        """Constructor."""
        super(PolygonizerDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://doc.qt.io/qt-5/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.printHelloWorld.clicked.connect(self.eventPushButtonDoSomethingOnClick)

    def compute_subsections(self, total_length, goal_length):
        if goal_length >= total_length:
            return 1, total_length
        else:
            num_subsections = round(total_length / goal_length)
            actual_length = total_length / num_subsections
            return num_subsections, actual_length

    def print_tx_geometry_as_geojson(self, geometry):
        source_crs = QgsCoordinateReferenceSystem("EPSG:2277")
        destination_crs = QgsCoordinateReferenceSystem("EPSG:4326")
        transform = QgsCoordinateTransform(source_crs, destination_crs, QgsProject.instance())
        geometry.transform(transform)
        print(geometry.asJson())

    def multiline_feature_to_linestring_geometry(self, feature):
        multiline = feature.geometry().asMultiPolyline()
        # Convert to LineString (it's up to you to not hand it a disjointed MultiLineString)
        points = [point for sublist in multiline for point in sublist]
        # Create a LineString from these points
        linestring = QgsLineString(points)
        return linestring

    def eventPushButtonDoSomethingOnClick(self):
        project = QgsProject.instance()

        layers = project.mapLayers()
        #print("Layers:")
        for layer_key in layers:
            layer = layers[layer_key]
            if layer.id().startswith("Workspace_"):
                #print("Layer", layer.id(), "is a", layer.type())
                project.removeMapLayer(layer.id())
        #print()

        # create layer
        selected_features_layer = QgsVectorLayer("LineString", "Workspace: Selected Features", "memory")
        crs = QgsCoordinateReferenceSystem("EPSG:2277")
        selected_features_layer.setCrs(crs)
        selected_features_provider = selected_features_layer.dataProvider()

        active_layer = iface.activeLayer()
        selected_features = active_layer.selectedFeatures()
        for intersection in selected_features:
            #print("Feature ID:", type(feature))
            #print("Feature ID:", feature.id())
            #print("Feature attributes:", feature.attributes())
            #print("Feature geometry:", feature.geometry().asWkt()) 
            selected_features_provider.addFeatures([intersection])
            #print()

        selected_features_layer.updateExtents()

        intersection_layer = QgsVectorLayer('Point?crs=EPSG:2277', 'Workspace: Intersection Points', 'memory')
        intersection_provider = intersection_layer.dataProvider()
        
        new_field = QgsField('legs', QVariant.String)
        intersection_provider.addAttributes([new_field])
        intersection_layer.updateFields()

        # iterate over features, find intersecting features, so we can create a point at each intersection
        features = selected_features_layer.getFeatures()
        list_of_features = list(features)
        for i in range(len(list_of_features)):
            for j in range(i+1, len(list_of_features)):
                #print(f"Comparing {list_of_features[i].id()} and {list_of_features[j].id()}")
                intersection = list_of_features[i].geometry().intersection(list_of_features[j].geometry())
                #print("Intersection:", intersection)
                intersection_feature = QgsFeature()
                intersection_feature.setGeometry(intersection)
                intersection_provider.addFeature(intersection_feature)
        
        intersection_layer.updateExtents()
        
        # Delete duplicate features, so we have a single intersection point for each intersection
        index = QgsSpatialIndex()
        delete_ids = []
        for intersection in intersection_layer.getFeatures():
            # If the feature's geometry is already in the index, then it's a duplicate
            if list(index.intersects(intersection.geometry().boundingBox())):
                # Store the feature's ID in our list of features to delete
                delete_ids.append(intersection.id())
            else:
                # If it's not in the index, add it now
                index.addFeature(intersection)
        #print("Delete IDs:", delete_ids)

        with edit(intersection_layer):
            intersection_layer.deleteFeatures(delete_ids)

        intersection_layer.updateExtents()

        segmented_roads_layer = QgsVectorLayer("LineString", "Workspace: Segments", "memory")
        crs = QgsCoordinateReferenceSystem("EPSG:2277")
        segmented_roads_layer.setCrs(crs)
        segmented_roads_provider = segmented_roads_layer.dataProvider()

        #leg_end_points_layer = QgsVectorLayer('Point?crs=EPSG:2277', 'Workspace: Leg Endpoints', 'memory')
        #leg_end_points_provider = leg_end_points_layer.dataProvider()

        intersection_polygons_layer = QgsVectorLayer("Polygon?crs=EPSG:2277", f"Workspace: Intersection Polygons", "memory")
        intersection_polygons_provider = intersection_polygons_layer.dataProvider()

        interconnect_polygons_layer = QgsVectorLayer("Polygon?crs=EPSG:2277", f"Workspace: Interconnect Polygons", "memory")
        interconnect_polygons_provider = intersection_polygons_layer.dataProvider()

        GOAL_SEGMENT_LENGTH = self.goalSegmentLength.value()
        GOAL_LEG_LENGTH = self.goalLegLength.value()
        BUFFER_LENGTH = self.polygonWidth.value()
        BUFFER_DETAIL = 20

        if False:
            for intersection in intersection_layer.getFeatures():
                print("Intersection:", intersection.id())
                for road in selected_features_layer.getFeatures():
                    if intersection.geometry().intersects(road.geometry()):
                        #print(f"Intersection {intersection.id()} intersects road {road.id()} of length {road.geometry().length()} feet.")

                        # this makes plain subsection features
                        subsections = self.compute_subsections(road.geometry().length(), GOAL_SEGMENT_LENGTH)
                        #print("Subsections: (count, practical_length)", subsections)
                        # this makes a feature layer of equally length segments.
                        for i in range(subsections[0]):
                            start_point = road.geometry().interpolate(i * subsections[1]).asPoint()
                            end_point = road.geometry().interpolate((i + 1) * subsections[1]).asPoint()
                            #print("Start point:", start_point)
                            #print("End point:", end_point)
                            start_distance = road.geometry().lineLocatePoint(QgsGeometry.fromPointXY(start_point))
                            end_distance = road.geometry().lineLocatePoint(QgsGeometry.fromPointXY(end_point))
                            #print("Start distance:", start_distance)
                            #print("End distance:", end_distance)
                            if start_distance > end_distance:
                                start_distance, end_distance = end_distance, start_distance
                            linestring = self.multiline_feature_to_linestring_geometry(road)
                            segment = linestring.curveSubstring(start_distance, end_distance)
                            segment_feature = QgsFeature()
                            segment_feature.setGeometry(QgsGeometry.fromPolyline(segment))
                            segmented_roads_provider.addFeatures([segment_feature])


        for intersection in intersection_layer.getFeatures():
            print()
            print("Intersection:", intersection.id())

            # create a layer for the leg segments
            leg_segments_layer = QgsVectorLayer("LineString", f"Workspace: Intersection {intersection.id()} Legs", "memory")
            crs = QgsCoordinateReferenceSystem("EPSG:2277")
            leg_segments_layer.setCrs(crs)
            leg_segments_provider = leg_segments_layer.dataProvider()

            leg_segments_buffer_layer = QgsVectorLayer("Polygon?crs=EPSG:2277", f"Workspace: Intersection {intersection.id()} Buffered Legs", "memory")
            leg_segments_buffer_provider = leg_segments_buffer_layer.dataProvider()

            legs = {}

            for road in selected_features_layer.getFeatures():
                if intersection.geometry().intersects(road.geometry()):
                    #print(f"Intersection {intersection.id()} intersects road {road.id()} of length {road.geometry().length()} feet.")
                    LEG_LENGTH = GOAL_LEG_LENGTH
                    # we have a road that is a leg, let's compute the length of the leg
                    # FIXME we need some snapping flexibility, because this can leave slivers
                    # FIXME we need to detect cul-de-sacs and extend the leg the length of them, possibly splitting it up
                    if road.geometry().length() < LEG_LENGTH * 2:
                        LEG_LENGTH = road.geometry().length() / 2

                    #print("Leg length:", LEG_LENGTH)

                    if False:
                        left_leg_end_point = road.geometry().interpolate(LEG_LENGTH).asPoint()
                        right_leg_end_point = road.geometry().interpolate(road.geometry().length() - LEG_LENGTH).asPoint()

                        left_leg_end_point_feature = QgsFeature()
                        left_leg_end_point_feature.setGeometry(QgsGeometry.fromPointXY(left_leg_end_point))
                        leg_end_points_provider.addFeature(left_leg_end_point_feature)

                        right_leg_end_point_feature = QgsFeature()
                        right_leg_end_point_feature.setGeometry(QgsGeometry.fromPointXY(right_leg_end_point))
                        leg_end_points_provider.addFeature(right_leg_end_point_feature)

                    linestring = self.multiline_feature_to_linestring_geometry(road)

                    left_leg  = linestring.curveSubstring(0, LEG_LENGTH)
                    right_leg = linestring.curveSubstring(linestring.length() - LEG_LENGTH, linestring.length())

                    left_leg_feature = QgsFeature()
                    left_leg_feature.setGeometry(QgsGeometry.fromPolyline(left_leg))

                    right_leg_feature = QgsFeature()
                    right_leg_feature.setGeometry(QgsGeometry.fromPolyline(right_leg))

                    if intersection.geometry().intersects(left_leg_feature.geometry()):
                        #print("Intersection intersects left leg")
                        leg_segments_provider.addFeatures([left_leg_feature])
                        legs[road.id()] = left_leg_feature.geometry().length()
                    if intersection.geometry().intersects(right_leg_feature.geometry()):
                        #print("Intersection intersects right leg")
                        leg_segments_provider.addFeatures([right_leg_feature])
                        legs[road.id()] = right_leg_feature.geometry().length()
                    
            # Start editing
            intersection_layer.startEditing()

            print("Intersection: ", intersection)
            # Set the feature's attribute
            intersection['legs'] = json.dumps(legs)
            print("Legs: ", intersection['legs'])

            # Update the feature in the layer
            #intersection_provider.updateFeature(intersection)
            intersection_layer.changeAttributeValue(intersection.id(), 0, intersection['legs'])

            #intersection.

            # Commit changes
            intersection_layer.commitChanges()
                    
            for leg in leg_segments_layer.getFeatures():
                #print(leg)

                end_cap_style = Qgis.EndCapStyle(2) # flat
                join_style = Qgis.JoinStyle(1) # miter
                miter_limit = 2

                buffer = leg.geometry().buffer(distance=BUFFER_LENGTH, segments=BUFFER_DETAIL, endCapStyle=end_cap_style, joinStyle=join_style, miterLimit=miter_limit)

                #self.print_tx_geometry_as_geojson(buffer)

                buffered_leg = QgsFeature()
                buffered_leg.setGeometry(buffer)
                leg_segments_buffer_provider.addFeature(buffered_leg)

            intersection_buffer = intersection.geometry().buffer(distance=BUFFER_LENGTH, segments=BUFFER_DETAIL)
            buffered_intersection = QgsFeature()
            buffered_intersection.setGeometry(intersection_buffer)
            leg_segments_buffer_provider.addFeature(buffered_intersection)
            #self.print_tx_geometry_as_geojson(intersection_buffer)


            features = leg_segments_buffer_layer.getFeatures()
            geometries = [feature.geometry() for feature in features]
            union_geometry = QgsGeometry.unaryUnion(geometries)
            union_feature = QgsFeature()
            union_feature.setGeometry(union_geometry)
            intersection_polygons_provider.addFeature(union_feature)

            leg_segments_layer.updateExtents()
            #project.addMapLayer(leg_segments_layer)

            leg_segments_buffer_layer.updateExtents()
            #project.addMapLayer(leg_segments_buffer_layer)

            intersection_polygons_layer.updateExtents()
            project.addMapLayer(intersection_polygons_layer)

        #leg_end_points_layer.updateExtents()
        #leg_segments_layer.updateExtents()
        
        project.addMapLayer(selected_features_layer)
        project.addMapLayer(intersection_layer)
        #project.addMapLayer(segmented_roads_layer)
        #project.addMapLayer(leg_end_points_layer)

        for road in selected_features_layer.getFeatures():
            print()
            print("Road:", road.id(), "Length:", road.geometry().length())

            linestring = self.multiline_feature_to_linestring_geometry(road)
            start_point = QgsGeometry.fromWkt(linestring.startPoint().asWkt())
            end_point = QgsGeometry.fromWkt(linestring.endPoint().asWkt())
            
            start_point_intersection_id = None
            end_point_intersection_id = None
            start_point_leg_length = None
            end_point_leg_length = None

            for intersection in intersection_layer.getFeatures():
                intersection_legs = json.loads(intersection["legs"])
                # print(f"Intersection {intersection.id()} legs:", intersection_legs)
                intersection_geometry = intersection.geometry()
                #print("Intersection Geometry:", intersection_geometry)
                #print("Start point:", start_point)
                if intersection_geometry.intersects(start_point):
                    start_point_intersection_id = intersection.id()
                    start_point_leg_length = intersection_legs[str(road.id())]
                if intersection_geometry.intersects(end_point):
                    end_point_intersection_id = intersection.id()
                    end_point_leg_length = intersection_legs[str(road.id())]

            start_point_leg_length = 0.0 if start_point_leg_length is None else start_point_leg_length
            end_point_leg_length = 0.0 if end_point_leg_length is None else end_point_leg_length

            print("Start point intersection id:", start_point_intersection_id)
            print("End point intersection id:", end_point_intersection_id)
            print("Start point leg length:", start_point_leg_length)
            print("End point leg length:", end_point_leg_length)

            interconnect_length = road.geometry().length() - start_point_leg_length - end_point_leg_length
            print("Interconnect length:", interconnect_length)
            if interconnect_length > 0:
                subsections = self.compute_subsections(interconnect_length, GOAL_SEGMENT_LENGTH)
                print("Subsections:", subsections)
            

                for i in range(subsections[0]):
                    start_point = road.geometry().interpolate(start_point_leg_length + (i * subsections[1])).asPoint()
                    end_point = road.geometry().interpolate(start_point_leg_length + ((i + 1) * subsections[1])).asPoint()
                    #print("Start point:", start_point)
                    #print("End point:", end_point)
                    start_distance = road.geometry().lineLocatePoint(QgsGeometry.fromPointXY(start_point))
                    end_distance = road.geometry().lineLocatePoint(QgsGeometry.fromPointXY(end_point))
                    #print("Start distance:", start_distance)
                    #print("End distance:", end_distance)
                    if start_distance > end_distance:
                        start_distance, end_distance = end_distance, start_distance
                    linestring = self.multiline_feature_to_linestring_geometry(road)
                    segment = linestring.curveSubstring(start_distance, end_distance)
                    segment_feature = QgsFeature()
                    segment_feature.setGeometry(QgsGeometry.fromPolyline(segment))
                    segmented_roads_provider.addFeatures([segment_feature])

                    buffer = segment_feature.geometry().buffer(distance=BUFFER_LENGTH, segments=BUFFER_DETAIL, endCapStyle=end_cap_style, joinStyle=join_style, miterLimit=miter_limit)
                    buffered_feature = QgsFeature()
                    buffered_feature.setGeometry(buffer)
                    interconnect_polygons_provider.addFeature(buffered_feature)


        segmented_roads_layer.updateExtents()
        project.addMapLayer(segmented_roads_layer)
        interconnect_polygons_layer.updateExtents()
        project.addMapLayer(interconnect_polygons_layer)
            
            # start_point_is_contained = False
            # end_point_is_contained = False
            # for intersection in intersection_layer.getFeatures():
                # target_point_geometry = intersection.geometry()
                # print("Intersection:", intersection["legs"])
                #print("Road geometry:",  road.geometry())

                
                #print("Start point:", start_point)
                #print("End point:", end_point)
                
                # if target_point_geometry.intersects(QgsGeometry.fromWkt(start_point.asWkt())):
                    # start_point_is_contained = True
                # if target_point_geometry.intersects(QgsGeometry.fromWkt(end_point.asWkt())):
                    # end_point_is_contained = True

            # print("Start point is contained:", start_point_is_contained)
            # print("End point is contained:", end_point_is_contained)

            # if start_point_is_contained and end_point_is_contained:
                # pass
            # elif start_point_is_contained and not end_point_is_contained:
                # pass
            # elif not start_point_is_contained and end_point_is_contained:
                # pass
            # else: # neither start nor end point is contained 🔥
                # pass

            # check each end for an intersection, if found 1, cul de sac, if found 2 an interconnect


    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()

