# -*- coding: utf-8 -*-
"""
/***************************************************************************
 *                                                                         *
 *  This software is released under the Unlicense. You may use it for      *
 *  anything you wish. This software is provided "AS-IS", and without      *
 *  warranty of any kind, express or implied.                              *
 *                                                                         *
 *  Please refer to the LICENSE file for more information.                 *
 *                                                                         *
 *  Additionally, please refer to <https://unlicense.org>                  *
 *                                                                         *
 ***************************************************************************/
"""

import os
import json

from qgis.PyQt import QtGui, QtWidgets, uic
from qgis.PyQt.QtCore import pyqtSignal, QVariant

from qgis.core import *
from qgis.utils import *

FORM_CLASS, _ = uic.loadUiType(
    os.path.join(os.path.dirname(__file__), "polygonizer_dockwidget_base.ui")
)


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
        self.doPolygonizerButton.clicked.connect(self.eventPushButtonRunPolygonizerOnClick)

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
        transform = QgsCoordinateTransform(
            source_crs, destination_crs, QgsProject.instance()
        )
        geometry.transform(transform)
        print(geometry.asJson())

    def multiline_feature_to_linestring_geometry(self, feature):
        multiline = feature.geometry().asMultiPolyline()
        # Convert to LineString (it's up to you to not hand it a disjointed MultiLineString)
        points = [point for sublist in multiline for point in sublist]
        # Create a LineString from these points
        linestring = QgsLineString(points)
        return linestring

    def add_layer_to_output_group(
        self, added_layer, output_group_name, layer_root, visible=True
    ):
        project = QgsProject.instance()
        layer_to_move = layer_root.findLayer(added_layer.id())
        if layer_to_move is not None:
            # Clone the layer
            clone = layer_to_move.clone()

            # Find the group
            group = layer_root.findGroup(output_group_name)

            # If the group was found, add the layer to it
            if group is not None:
                # Add the cloned layer to the group
                group.insertChildNode(-1, clone)
                if visible:
                    clone.setItemVisibilityChecked(True)
                else:
                    clone.setItemVisibilityChecked(False)

                # Remove the original layer from its parent
                parent = layer_to_move.parent()
                parent.removeChildNode(layer_to_move)

    def eventPushButtonRunPolygonizerOnClick(self):
        project = QgsProject.instance()

        # clean up previous executions workspace layers
        layers = project.mapLayers()
        for layer_key in layers:
            layer = layers[layer_key]
            if layer.id().startswith("Workspace_"):
                project.removeMapLayer(layer.id())

        # create a group and move it to the top of the layer list for output
        output_group_name = "Polygonizer Output"
        layer_root = project.layerTreeRoot()
        group = layer_root.findGroup(output_group_name)
        # If the group was found, remove it
        if group is not None:
            parent = group.parent()
            parent.removeChildNode(group)
        output_group = layer_root.addGroup(output_group_name)
        clone = output_group.clone()
        layer_root.insertChildNode(0, clone)
        parent = output_group.parent()
        parent.removeChildNode(output_group)

        # create a layer of the selected features
        selected_features_layer = QgsVectorLayer(
            "LineString", "Workspace: Selected Features", "memory"
        )
        crs = QgsCoordinateReferenceSystem("EPSG:2277")
        selected_features_layer.setCrs(crs)
        selected_features_provider = selected_features_layer.dataProvider()
        active_layer = iface.activeLayer()
        selected_features = active_layer.selectedFeatures()
        for intersection in selected_features:
            selected_features_provider.addFeature(intersection)
        selected_features_layer.updateExtents()

        # create a layer of the intersection points of road segments, meaning intersections
        intersection_layer = QgsVectorLayer(
            "Point?crs=EPSG:2277", "Workspace: Intersection Points", "memory"
        )
        intersection_provider = intersection_layer.dataProvider()
        new_field = QgsField("legs", QVariant.String)
        intersection_provider.addAttributes([new_field])
        intersection_layer.updateFields()

        # iterate over features, find intersecting features, so we can create a point at each intersection
        features = selected_features_layer.getFeatures()
        list_of_features = list(features)
        for i in range(len(list_of_features)):
            for j in range(i + 1, len(list_of_features)):
                intersection = (
                    list_of_features[i]
                    .geometry()
                    .intersection(list_of_features[j].geometry())
                )
                intersection_feature = QgsFeature()
                intersection_feature.setGeometry(intersection)
                intersection_provider.addFeature(intersection_feature)

        # housekeeping on the layer we just populated
        intersection_layer.updateExtents()

        # find duplicated features and keep track of which we can delete
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

        # delete the duplicate features
        with edit(intersection_layer):
            intersection_layer.deleteFeatures(delete_ids)

        # housekeeping again
        intersection_layer.updateExtents()

        # create a layer to hold road segments between legs of an intersection
        segmented_roads_layer = QgsVectorLayer(
            "LineString", "Workspace: Segments", "memory"
        )
        crs = QgsCoordinateReferenceSystem("EPSG:2277")
        segmented_roads_layer.setCrs(crs)
        segmented_roads_provider = segmented_roads_layer.dataProvider()

        # create a layer of at the points where legs end from an intersection
        # leg_end_points_layer = QgsVectorLayer('Point?crs=EPSG:2277', 'Workspace: Leg Endpoints', 'memory')
        # leg_end_points_provider = leg_end_points_layer.dataProvider()

        # a layer to hold polygons we create for intersections
        intersection_polygons_layer = QgsVectorLayer(
            "Polygon?crs=EPSG:2277", f"Workspace: Intersection Polygons", "memory"
        )
        intersection_polygons_provider = intersection_polygons_layer.dataProvider()

        # a layer to hold polygons we create for interconnects between intersections
        interconnect_polygons_layer = QgsVectorLayer(
            "Polygon?crs=EPSG:2277", f"Workspace: Interconnect Polygons", "memory"
        )
        interconnect_polygons_provider = interconnect_polygons_layer.dataProvider()

        # ingest parameters of the polygonization algorithm
        GOAL_SEGMENT_LENGTH = self.idealSegmentLengthSpinbox.value()
        GOAL_LEG_LENGTH = self.idealLegLengthSpinbox.value()
        BUFFER_LENGTH = self.polygonWidthSpinbox.value()
        MAX_SNAP_LENGTH = self.maxSnapLengthSpinbox.value()
        BUFFER_DETAIL = 20

        end_cap_style = Qgis.EndCapStyle(2)  # flat
        join_style = Qgis.JoinStyle(1)  # miter
        miter_limit = 2

        # begin polygon algorithm and begin by iterating over the intersection center points
        for intersection in intersection_layer.getFeatures():

            # create a layer for the leg segments linestrings
            leg_segments_layer = QgsVectorLayer(
                "LineString",
                f"Workspace: Intersection {intersection.id()} Legs",
                "memory",
            )
            crs = QgsCoordinateReferenceSystem("EPSG:2277")
            leg_segments_layer.setCrs(crs)
            leg_segments_provider = leg_segments_layer.dataProvider()

            # create a layer for the buffered leg segment polygons that we'll union into the intersection polygon
            # we're making one of these layers per intersection
            leg_segments_buffer_layer = QgsVectorLayer(
                "Polygon?crs=EPSG:2277",
                f"Workspace: Intersection {intersection.id()} Buffered Legs",
                "memory",
            )
            leg_segments_buffer_provider = leg_segments_buffer_layer.dataProvider()

            legs = {}

            # iterate over roads, looking for ones which intersect with the polygon center points
            for road in selected_features_layer.getFeatures():
                if intersection.geometry().intersects(road.geometry()):
                    LEG_LENGTH = GOAL_LEG_LENGTH
                    # we have a road that is a leg, let's compute the length of the leg
                    if road.geometry().length() < LEG_LENGTH * 2:
                        LEG_LENGTH = road.geometry().length() / 2
                    elif (road.geometry().length() - (LEG_LENGTH * 2)) < MAX_SNAP_LENGTH:
                        LEG_LENGTH = road.geometry().length() / 2

                    # we don't know the start-to-stop orientation of the line, so let's grab some endpoints and look for intersections
                    linestring = self.multiline_feature_to_linestring_geometry(road)

                    left_leg = linestring.curveSubstring(0, LEG_LENGTH)
                    right_leg = linestring.curveSubstring(
                        linestring.length() - LEG_LENGTH, linestring.length()
                    )

                    left_leg_feature = QgsFeature()
                    left_leg_feature.setGeometry(QgsGeometry.fromPolyline(left_leg))

                    right_leg_feature = QgsFeature()
                    right_leg_feature.setGeometry(QgsGeometry.fromPolyline(right_leg))

                    # here we're going to do almost the same thing but different for whichever end of the line is the one in question
                    if intersection.geometry().intersects(left_leg_feature.geometry()):
                        leg_segments_provider.addFeatures([left_leg_feature])
                        legs[road.id()] = left_leg_feature.geometry().length()
                    if intersection.geometry().intersects(right_leg_feature.geometry()):
                        leg_segments_provider.addFeatures([right_leg_feature])
                        legs[road.id()] = right_leg_feature.geometry().length()

            # Start editing
            intersection_layer.startEditing()

            # turns out, we need to keep track of some metadata about leg lengths per road emanating from an intersection
            intersection["legs"] = json.dumps(legs)

            # we'll store that metadata in the intersection layer
            intersection_layer.changeAttributeValue(
                intersection.id(), 0, intersection["legs"]
            )

            # Commit changes
            intersection_layer.commitChanges()

            # iterate over the legs of the intersection and buffer them
            for leg in leg_segments_layer.getFeatures():
                buffer = leg.geometry().buffer(
                    distance=BUFFER_LENGTH,
                    segments=BUFFER_DETAIL,
                    endCapStyle=end_cap_style,
                    joinStyle=join_style,
                    miterLimit=miter_limit,
                )
                buffered_leg = QgsFeature()
                buffered_leg.setGeometry(buffer)
                leg_segments_buffer_provider.addFeature(buffered_leg)

            # we also need a buffer around the center point to fill in any slivers that the legs do not cover
            intersection_buffer = intersection.geometry().buffer(
                distance=BUFFER_LENGTH, segments=BUFFER_DETAIL
            )
            buffered_intersection = QgsFeature()
            buffered_intersection.setGeometry(intersection_buffer)
            leg_segments_buffer_provider.addFeature(buffered_intersection)
            # self.print_tx_geometry_as_geojson(intersection_buffer)

            # pull the features together that we've accumulated in this intersections buffer layer
            features = leg_segments_buffer_layer.getFeatures()
            geometries = [feature.geometry() for feature in features]
            union_geometry = QgsGeometry.unaryUnion(geometries)
            union_feature = QgsFeature()
            union_feature.setGeometry(union_geometry)
            intersection_polygons_provider.addFeature(union_feature)

            # housekeeping on the modified layers
            leg_segments_layer.updateExtents()
            leg_segments_buffer_layer.updateExtents()

        # housekeeping on some layers
        intersection_polygons_layer.updateExtents()
        selected_features_layer.updateExtents()
        intersection_layer.updateExtents()

        # add them to the map
        added_layer = project.addMapLayer(intersection_polygons_layer)
        self.add_layer_to_output_group(
            added_layer, output_group_name, layer_root, visible=True
        )

        added_layer = project.addMapLayer(selected_features_layer)
        self.add_layer_to_output_group(
            added_layer, output_group_name, layer_root, visible=False
        )

        added_layer = project.addMapLayer(intersection_layer)
        self.add_layer_to_output_group(
            added_layer, output_group_name, layer_root, visible=False
        )

        # back out here in the main function execution path, we're going to build the polygons for the intermediate sections between intersections
        for road in selected_features_layer.getFeatures():
            # get some raw geometry to grab endpoints from
            linestring = self.multiline_feature_to_linestring_geometry(road)
            start_point = QgsGeometry.fromWkt(linestring.startPoint().asWkt())
            end_point = QgsGeometry.fromWkt(linestring.endPoint().asWkt())

            # places to keep data
            start_point_intersection_id = None
            end_point_intersection_id = None
            start_point_leg_length = None
            end_point_leg_length = None

            # let's run over the intersections, looking for ones that are involved in this road
            for intersection in intersection_layer.getFeatures():
                intersection_legs = json.loads(intersection["legs"])
                intersection_geometry = intersection.geometry()
                # found one, let's grab the intersection leg length
                if intersection_geometry.intersects(start_point):
                    start_point_intersection_id = intersection.id()
                    start_point_leg_length = intersection_legs[str(road.id())]
                # and found the other one
                if intersection_geometry.intersects(end_point):
                    end_point_intersection_id = intersection.id()
                    end_point_leg_length = intersection_legs[str(road.id())]

            # a road may end at a dead end, so we need to handle that case
            start_point_leg_length = (
                0.0 if start_point_leg_length is None else start_point_leg_length
            )
            end_point_leg_length = (
                0.0 if end_point_leg_length is None else end_point_leg_length
            )

            # finally compute the total length of the non-intersectional interconnect portion of the road
            interconnect_length = (
                road.geometry().length() - start_point_leg_length - end_point_leg_length
            )
            # GIS & floats: the age old problem. Define a cut off for what we believe is zero.
            # This is about 2500 hydrogen atoms long. This is very zero in the scope of a road.
            if (
                interconnect_length > 0.0000001
            ):  
                # calculate the number of subsections we'll need to create
                subsections = self.compute_subsections(
                    interconnect_length, GOAL_SEGMENT_LENGTH
                )

                # get the geometry we're going to be dicing up
                road_geometry = road.geometry()

                # iterate over every subsection of intersectional road we want to grab a segment of to buffer
                for i in range(subsections[0]):
                    start_point = (
                        road.geometry()
                        .interpolate(start_point_leg_length + (i * subsections[1]))
                        .asPoint()
                    )
                    try:
                        end_point = (
                            road.geometry()
                            .interpolate(
                                start_point_leg_length + ((i + 1) * subsections[1])
                            )
                            .asPoint()
                        )
                    # handle this edge case where zero didn't register as zero and you have a ultra tiny segment
                    except:
                        end_point = (
                            road.geometry()
                            .interpolate(
                                start_point_leg_length
                                + ((i + 1) * subsections[1])
                                - 0.0000001
                            )
                            .asPoint()
                        )

                    # get some distances out of the road length and our determined point along the line
                    start_distance = road.geometry().lineLocatePoint(
                        QgsGeometry.fromPointXY(start_point)
                    )
                    end_distance = road.geometry().lineLocatePoint(
                        QgsGeometry.fromPointXY(end_point)
                    )
                    # make sure our math is sane
                    if start_distance > end_distance:
                        start_distance, end_distance = end_distance, start_distance
                    # get raw geometry to carve out the subsection of road linestring
                    linestring = self.multiline_feature_to_linestring_geometry(road)
                    # this line extracts the particular subsegment of the linestring we're interested in
                    segment = linestring.curveSubstring(start_distance, end_distance)
                    # and store that feature that was created
                    segment_feature = QgsFeature()
                    segment_feature.setGeometry(QgsGeometry.fromPolyline(segment))
                    segmented_roads_provider.addFeatures([segment_feature])

                    # and buffer that feature
                    buffer = segment_feature.geometry().buffer(
                        distance=BUFFER_LENGTH,
                        segments=BUFFER_DETAIL,
                        endCapStyle=end_cap_style,
                        joinStyle=join_style,
                        miterLimit=miter_limit,
                    )

                    # these are two cases to handle when the road ends not at an intersection, one for each potential non-intersectional endpoint
                    # we're going to buffer that point and join it to the last polygon so we get a faked round end cap but only on that one end
                    if not start_point_intersection_id and i == 0:
                        buffered_end = QgsGeometry.fromWkt(start_point.asWkt()).buffer(
                            distance=BUFFER_LENGTH, segments=BUFFER_DETAIL
                        )
                        buffer = buffer.combine(buffered_end)
                    if not end_point_intersection_id and i == subsections[0] - 1:
                        buffered_end = QgsGeometry.fromWkt(end_point.asWkt()).buffer(
                            distance=BUFFER_LENGTH, segments=BUFFER_DETAIL
                        )
                        buffer = buffer.combine(buffered_end)

                    # store the buffered feature in a layer of intersectional polygons
                    buffered_feature = QgsFeature()
                    buffered_feature.setGeometry(buffer)
                    interconnect_polygons_provider.addFeature(buffered_feature)

        # housekeeping for our edited layers
        segmented_roads_layer.updateExtents()
        interconnect_polygons_layer.updateExtents()

        # add our new layer with our intersectional polygons to the output group
        interconnect_added_layer = project.addMapLayer(interconnect_polygons_layer)
        self.add_layer_to_output_group(
            interconnect_added_layer, output_group_name, layer_root
        )

    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()
