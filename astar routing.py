import sys
import numpy as np
import math
# -----------------------------
import pandas as pd
# import matplotlib.pyplot as plt

# ------------------------------
#import csv
import networkx as nx

from qgis.core import QgsApplication, QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsPoint, \
    QgsProcessingFeedback
from qgis.analysis import QgsNativeAlgorithms

QgsApplication.setPrefixPath(sys.prefix, True)
sys.path.append('/Program Files/QGIS 3.34.6/apps/qgis-ltr/python/plugins')
from processing.core.Processing import Processing
import processing

#import time

# -------------------------------
qgs = QgsApplication([], True)
qgs.initQgis()
qgs.processingRegistry().addProvider(QgsNativeAlgorithms())
Processing.initialize()

registry = QgsProject.instance()
registry.read('/design_optim.qgz')
points = registry.mapLayersByName('grid_slope')
feats_points = [feat for feat in points[0].getFeatures()]
overlay = registry.mapLayersByName('area_geo')[0]
# ----------------cost-----------------
stringing_cost = {26: 52, 36: 83}  # $/m
welding_cost = {26: 3204, 36: 10664}  # $/Num
welding_inspection_cost = {26: 1180, 36: 1573}  # $/Num
coating_cost = {26: 251, 36: 440}  # $/m
rbore_construction_cost = 285621  # $/Num
rcross_construction_cost = {26: 173490, 36: 238067}  # $/Num
hydrotesting_cost = {26: 86, 36: 129}  # $/m
misc_hydrotesting_cost = {26: 7327, 36: 11235}  # $/km
cat_protection_testing_cost = 2164  # $/km
cat_protection_installation_cost = 2856  # $/Num
pipe_hauling_cost = 165  # $/ton
trenching_cost = {'light_soil': 18, 'rocky_surface': 199}  # $/m^3
backfilling_cost = 13  # $/m^3
comp_install_cost = 2500000  # $  (sanaye 2013)
cf_river = 1  # to be calculated
cf_road = 1
n_weight = {26: 8, 36: 3}
trenching_width = {26: 32, 36: 44}
trenching_depth = 2
# -------------------------------------
numx = 170
numy = 141

G = nx.grid_2d_graph(numx, numy)

G.add_edges_from([
                     ((x, y), (x + 1, y + 1))
                     for x in range(numx - 1)
                     for y in range(numy - 1)
                 ] + [
                     ((x + 1, y), (x, y + 1))
                     for x in range(numx - 1)
                     for y in range(numy - 1)
                 ])
pos = {(x, y): (x, -y) for x, y in G.nodes()}
with pd.ExcelFile('/grid_points.xlsx') as xls:
    df1 = pd.read_excel(xls, 'RemovedPoints')
    for elem in df1['id']:
        x = (elem // numy) * (elem % numy != 0) + (elem // numy - 1) * (elem % numy == 0)
        y = (elem % numy - 1) * (elem % numy != 0) + (numy - 1) * (elem % numy == 0)
        G.remove_node((x, y))
# ----------------problem input
diameter = {1: 26, 2: 36, 3: 26, 4: 36}


# -----------------------------
class AstarNetworkRouting:
    def __init__(self):
        self.expansion_cache = {}
        self.diam = 0
        self.distance_based = 0

    def reset_assign(self):
        pass

    def _in_exp_cache(self, diam, start_node, end_node):
        return False

    def dist(self, a, b):
        (x1, y1) = a
        (x2, y2) = b
        pid1 = x1 * numy + y1
        pid2 = x2 * numy + y2
        PointFeatList = []
        PointFeatList.append(feats_points[pid1])
        PointFeatList.append(feats_points[pid2])

        PointList = []
        for elem in PointFeatList:
            temp_point = QgsPoint(float(elem.attribute('left')), float(elem.attribute('top')))
            PointList.append(temp_point)

        distance = PointList[0].distance(PointList[1])
        # distance = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
        estimated_cost = 61.1 - 5.75 * 0.000621371 * distance - 1.09 * self.diam + \
                         0.36 * 0.000621371 * distance * self.diam

        if self.distance_based == 1:
            return distance
        else:
            return estimated_cost * 1e6

    def cost_func(self, p1, p2, e_att):
        (x1, y1) = p1
        (x2, y2) = p2
        pid1 = x1 * numy + y1
        pid2 = x2 * numy + y2
        # -------------calculate slope and cost factor
        s1 = feats_points[pid1].attribute('slope')
        s2 = feats_points[pid2].attribute('slope')
        slope = (s1 + s2) / 2
        cf_general = 1 * (slope <= 8) + 1.15 * (8 <= slope < 15) + 1.2 * (15 <= slope < 22) + 1.4 * (slope >= 22)
        cf_hydrotest = 1 * (slope <= 8) + 1.25 * (8 <= slope < 15) + 1.30 * (15 <= slope < 22) + 1.45 * (
                slope >= 22)
        # ------------calculate rocky surface passage
        PointFeatList = []
        PointFeatList.append(feats_points[pid1])
        PointFeatList.append(feats_points[pid2])

        PointList = []
        for elem in PointFeatList:
            temp_point = QgsPoint(float(elem.attribute('left')), float(elem.attribute('top')))
            PointList.append(temp_point)

        vl = QgsVectorLayer(
            "LineString?crs=epsg:3857&field=id:integer&index=yes",
            'candidateline', 'memory')

        vl.startEditing()
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry.fromPolyline(PointList))
        feature.setAttributes([1])
        vl.addFeature(feature)
        vl.commitChanges()
        route_length = 0.0
        for f in vl.getFeatures():
            route_length = f.geometry().length()
        params = {'INPUT': vl, 'OVERLAY': overlay, 'OUTPUT': 'TEMPORARY_OUTPUT'}
        feedback = QgsProcessingFeedback()
        res = processing.run("native:difference", params, feedback=feedback)
        dl = res['OUTPUT']
        rocky_length = route_length
        for f in dl.getFeatures():
            rocky_length = route_length - f.geometry().length()

        pipe_num = route_length // 13 + 1
        welding_point_num = pipe_num - 1
        total_pipe_weight = 28.8 * pipe_num / n_weight[self.diam]
        cat_protect_num = (route_length / 1000) // 10 + 1

        str_costs = stringing_cost[self.diam] * route_length * cf_general
        wel_costs = welding_cost[self.diam] * welding_point_num * cf_general
        wel_ins_costs = welding_inspection_cost[self.diam] * welding_point_num
        coat_costs = coating_cost[self.diam] * route_length * cf_general
        rb_costs = 0
        rc_costs = 0
        hydtest_costs = hydrotesting_cost[self.diam] * route_length * cf_hydrotest
        mischydro_costs = misc_hydrotesting_cost[self.diam] * route_length / 1000
        cptest_costs = cat_protection_testing_cost * route_length / 1000
        cpi_costs = cat_protection_installation_cost * cat_protect_num
        pipeh_costs = pipe_hauling_cost * total_pipe_weight
        tr_costs = trenching_width[self.diam] * 25.4 / 1000 * trenching_depth \
                   * (trenching_cost['light_soil'] * (route_length - rocky_length) + trenching_cost[
            'rocky_surface'] * rocky_length)
        bfill_costs = backfilling_cost * trenching_width[self.diam] * 25.4 / 1000 * trenching_depth * route_length

        total_costs = (str_costs + wel_costs + wel_ins_costs + coat_costs + rb_costs + rc_costs + hydtest_costs
                       + mischydro_costs + cptest_costs + cpi_costs + pipeh_costs + tr_costs + bfill_costs)
        mod_total_costs = total_costs / (1 - .39)

        return mod_total_costs

    def preparation(self, expansion_result):
        pass

    def operation(self, demand_factor):
        pass

    def expansion(self, diam, comp_enabled, start_node, end_node):
        if not self._in_exp_cache(diam, start_node, end_node):
            self.diam = diam  # A* cost func does not accept extra arguments
            cfunc = lambda p1, p2, e_att: self.cost_func(p1, p2, e_att)
            hfunc = lambda a, b: self.dist(a, b)
            opt_path = nx.astar_path(G, start_node, end_node, heuristic=hfunc, weight=cfunc)
            opt_path_cost = sum(cfunc(u, v, G[u][v]) for u, v in zip(opt_path[:-1], opt_path[1:]))
            self.distance_based = 1
            opt_path_length = sum(hfunc(a, b) for a, b in zip(opt_path[:-1], opt_path[1:]))
            self.distance_based = 0
            self.expansion_cache[(diam, start_node, end_node)] = [opt_path, opt_path_length, opt_path_cost]
        opt_path, opt_path_length, opt_path_cost = self.expansion_cache[(diam, start_node, end_node)]
        opt_path_cost += comp_enabled * comp_install_cost  # add compressor installation costtwo
        return opt_path, opt_path_length, opt_path_cost

    def fill_array(self, a, n):
        pass

    def exp_cons(self, x):
        pass

    def execute(self, x):

        diam, comp_enabled, start_node, end_node = x
        opt_path, opt_path_length, opt_path_cost = self.expansion(diam, comp_enabled, start_node, end_node)
        expansion_result = [opt_path, opt_path_length, opt_path_cost]

        return expansion_result


exp_test = AstarNetworkRouting()
input = []
exp_res = exp_test.execute(input)
