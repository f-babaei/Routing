import numpy as np
import math
from scipy.optimize import NonlinearConstraint, Bounds, differential_evolution
# -----------------------------
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------
import csv
import networkx as nx

from qgis.core import QgsApplication, QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsPoint, \
    QgsProcessingFeedback
from qgis.analysis import QgsNativeAlgorithms

QgsApplication.setPrefixPath('/usr', True)

from processing.core.Processing import Processing
import processing
import time

# -------------------------------

qgs = QgsApplication([], True)
qgs.initQgis()
qgs.processingRegistry().addProvider(QgsNativeAlgorithms())
Processing.initialize()

registry = QgsProject.instance()
registry.read('/home/controllab/Downloads/design_optim.qgz')
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
with pd.ExcelFile('/home/controllab/Downloads/grid_points.xlsx') as xls:
    df1 = pd.read_excel(xls, 'RemovedPoints')
    for elem in df1['id']:
        x = (elem // numy) * (elem % numy != 0) + (elem // numy - 1) * (elem % numy == 0)
        y = (elem % numy - 1) * (elem % numy != 0) + (numy - 1) * (elem % numy == 0)
        G.remove_node((x, y))
# ----------------problem input
diameter = {1: 26, 2: 36, 3: 26, 4: 36}
nInit = 6
nInter = 5
nTarget = 3
nN = nInit + nInter + nTarget
'''with pd.ExcelFile('/home/controllab/Downloads/tech_flow.xlsx') as xls:
    tech_flow = pd.read_excel(xls, 'Sheet1')'''

opt_to_graph = {0: (41, 137), 1: (88, 73), 2: (63, 55), 3: (40, 77), 4: (169, 58), 5: (169, 9),
                6: (82, 110), 7: (24, 103), 8: (122, 80), 9: (113, 59), 10: (154, 44),
                11: (95, 118), 12: (0, 97), 13: (125, 48)}
node_neighbors = {0: [6, 7, 11, 12], 1: [6, 8, 9, 11, 13], 2: [9, 13], 3: [7, 12], 4: [8, 10, 11, 13], 5: [10, 13],
                  6: [8, 11], 7: [12], 8: [11], 9: [10, 13], 10: [13],
                  11: [], 12: [], 13: []}
node_priority = {11: (0, 1, 4), 12: (3, 0), 13: (2, 1, 4, 5)}
aux_map = {(0, 1): [6], (0, 1, 4): [6, 8], (3, 0): [7], (2, 1): [9], (2, 1, 4, 5): [9, 10]}
nT = 0
for elem in node_neighbors:
    nT += len(node_neighbors[elem])
nT = 12
bigM = 2

obj_list = []
conv_list = []


# -----------------------------
class my_minlp:
    def __init__(self):
        self.MATeng = me.start_matlab()
        self.counter = 0
        self.fevalcount = 0
        self.diam = 0
        self.distance_based = 0

        self.assignable_node_id = 34
        self.assignable_pipe_id = 25
        self.assignable_comp_id = 6
        self.installed_pipe_id = []
        self.expansion_cache = {}
        self.installed_node = [0, 1, 2, 3, 4, 5, 11, 12, 13]
        self.opt_to_grail = {0: 9, 1: 10, 2: 12, 3: 13, 4: 14, 5: 22, 6: None, 7: None, 8: None, 9: None, 10: None,
                             11: 31, 12: 32, 13: 33}

    def reset_assign(self):
        self.assignable_node_id = 34
        self.assignable_pipe_id = 25
        self.assignable_comp_id = 6
        self.installed_pipe_id = []
        self.installed_node = [0, 1, 2, 3, 4, 5, 11, 12, 13]
        self.opt_to_grail = {0: 9, 1: 10, 2: 12, 3: 13, 4: 14, 5: 22, 6: None, 7: None, 8: None, 9: None, 10: None,
                             11: 31, 12: 32, 13: 33}

    def _in_exp_cache(self, diam, start_node, end_node):
        return ((diam, start_node, end_node) in self.expansion_cache) or (
                (diam, end_node, start_node) in self.expansion_cache)

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
        estimated_cost = 61.1 - 5.75 * 0.000621371 * distance - 1.09 * self.diam + 0.36 * 0.000621371 * distance * self.diam
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
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/original_data/model30t_1'
                  '/input_network_comps.csv', newline='') as csv_file:
            comp_data = [line for line in csv.reader(csv_file)]
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/original_data/model30t_1'
                  '/input_network_pipes.csv', newline='') as csv_file:
            pipe_data = [line for line in csv.reader(csv_file)]
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/original_data/model30t_1'
                  '/input_network_nodes.csv', newline='') as csv_file:
            node_data = [line for line in csv.reader(csv_file)]
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/original_data/model30t_1'
                  '/input_ts_qbar.csv', newline='') as csv_file:
            qbar_data = [line for line in csv.reader(csv_file)]

        for elem in expansion_result:
            if not self.opt_to_grail[elem[0]]:
                self.opt_to_grail[elem[0]] = self.assignable_node_id
                self.assignable_node_id += 1
            if not self.opt_to_grail[elem[1]]:
                self.opt_to_grail[elem[1]] = self.assignable_node_id
                self.assignable_node_id += 1

            start_node_id = self.opt_to_grail[elem[0]]
            end_node_id = self.opt_to_grail[elem[1]]
            for item in elem:
                if item not in self.installed_node:
                    new_node = []
                    new_node.append(self.opt_to_grail[item])  # id
                    new_node.append('Node{}'.format(self.opt_to_grail[item]))  # name
                    new_node.append(0.8714)  # x coord
                    new_node.append(-.755)  # y coord
                    new_node.append(500)  # min p
                    new_node.append(800)  # max p
                    new_node.append(0)  # min injection
                    new_node.append(0)  # max injection
                    new_node.append(0)  # slack bool
                    node_data.append(new_node)
                    self.installed_node.append(item)

                    temp_z = [item[0] for item in np.zeros([102, 1], int).tolist()]
                    temp_z[0] = new_node[0]
                    for item in qbar_data:
                        item.append(temp_z[qbar_data.index(item)])

            operation_length = expansion_result[elem][0] / 1000 * 0.621371  # m to mile
            operation_diam = expansion_result[elem][1]
            operation_comp = expansion_result[elem][2]

            from_node_pipe = start_node_id

            if operation_comp == 1:
                new_node = []
                new_node.append(self.assignable_node_id)  # id
                new_node.append('Node{}'.format(self.assignable_node_id))  # name
                new_node.append(0.8714)  # x coord
                new_node.append(-.755)  # y coord
                new_node.append(500)  # min p
                new_node.append(800)  # max p
                new_node.append(0)  # min injection
                new_node.append(0)  # max injection
                new_node.append(0)  # slack bool
                self.assignable_node_id += 1
                node_data.append(new_node)

                temp_z = [item[0] for item in np.zeros([102, 1], int).tolist()]
                temp_z[0] = new_node[0]
                for item in qbar_data:
                    item.append(temp_z[qbar_data.index(item)])

                new_comp = []
                new_comp.append(self.assignable_comp_id)  # id
                new_comp.append('Comp{}'.format(self.assignable_comp_id))  # name
                new_comp.append(start_node_id)  # from node
                new_comp.append(new_node[0])  # to node
                new_comp.append(1)  # c min
                new_comp.append(1.6)  # c max
                new_comp.append(1500)  # power max # can be decision variable based on comp type
                new_comp.append(0)  # flow min
                new_comp.append(600)  # flow max
                from_node_pipe = new_comp[3]  # to node comp
                self.assignable_comp_id += 1
                comp_data.append(new_comp)

            new_pipe = []
            new_pipe.append(self.assignable_pipe_id)  # id
            self.installed_pipe_id.append(self.assignable_pipe_id)
            new_pipe.append('Pipe{}'.format(self.assignable_pipe_id))  # name
            new_pipe.append(from_node_pipe)  # from node
            new_pipe.append(end_node_id)  # to node
            new_pipe.append(operation_diam)  # diameter
            new_pipe.append(operation_length)  # length
            new_pipe.append(.01)  # friction
            new_pipe.append(0)  # disc seg
            self.assignable_pipe_id += 1
            pipe_data.append(new_pipe)

        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/data/model30t_1'
                  '/input_network_pipes.csv', 'w', newline='') as csv_file:
            csv.writer(csv_file).writerows(pipe_data)
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/data/model30t_1'
                  '/input_network_comps.csv', 'w', newline='') as csv_file:
            csv.writer(csv_file).writerows(comp_data)
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/data/model30t_1'
                  '/input_network_nodes.csv', 'w', newline='') as csv_file:
            csv.writer(csv_file).writerows(node_data)
        with open('/home/controllab/PycharmProjects/pythonProject/grail-master-design_two/data/model30t_1'
                  '/input_ts_qbar.csv', 'w', newline='') as csv_file:
            csv.writer(csv_file).writerows(qbar_data)

    def operation(self, demand_factor):

        gas_op_costs = self.MATeng.optim_grail(demand_factor)

        return gas_op_costs

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
        mask = np.tri(n, dtype=bool, k=-1)  # or np.arange(n)[:,None] > np.arange(n)
        mask[:, :] = False
        for elem in node_neighbors:
            for item in node_neighbors[elem]:
                mask[elem, item] = True
        out = np.zeros((n, n), dtype=int)
        out[mask] = a
        return out

    def exp_cons(self, x):

        yp = np.array([x[0], x[4], x[8]])
        ya = np.array([x[1], x[5], x[9]])
        yd = np.array([x[2], x[6], x[10]])
        yc = np.array([x[3], x[7], x[11]])

        c = np.sum(yp) - 1  # there exists a connection between init-intermediate and target
        c = np.append(c, yp - np.multiply(ya, bigM))
        c = np.append(c, yd - yp)
        c = np.append(c, yp - yc)
        c = np.append(c, np.multiply(yp, bigM) - yd)
        c = np.append(c, (yp[-1] != 3) - 1)

        return c.T

    '''def callback_de(self, xk, val):
        if self.exp_cons(xk)>=0:
            obj_list.append(self.fitness(xk))
        else:
            obj_list.append(10**11)
		
        conv_list.append(val)	
        return false'''

    def fitness(self, x):
        installed_comp = [0, 4]
        expansion_costs = 0
        expansion_result = {}

        '''[n-3:n] technology of valorization based on three demands (small, nominal, large)
        with open('/home/farshid/PycharmProjects/pythonProject/grail-master-design/original_data/model30t_1'
                  '/input_ts_gbar.csv', newline='') as csv_file:
            gbar_data = [line for line in csv.reader(csv_file)]
            
              for elem in f:
                  temp_z = tech_flow[item].tolist()
                for item in gbar_data:
                    item.append(temp_z[gbar_data.index(item)])
            
        with open('/home/farshid/PycharmProjects/pythonProject/grail-master-design/data/model30t_1'
          '/input_ts_gbar.csv', 'w', newline='') as csv_file:
          csv.writer(csv_file).writerows(gbar_data)
        '''

        yp = np.array([x[0], x[4], x[8]], dtype=int)
        ya = np.array([x[1], x[5], x[9]], dtype=int)
        yd = np.array([x[2], x[6], x[10]], dtype=int)
        yc = np.array([x[3], x[7], x[11]], dtype=int)
        for i in range(yp.size):
            start_nodes = node_priority[i + 11][0:yp[i]]
            if ya[i] == 1:
                aux_points = aux_map[start_nodes]
                for item in aux_points:
                    for elem in start_nodes:
                        if (item in node_neighbors[elem]) or (elem in node_neighbors[item]):
                            start_node = opt_to_graph[elem]
                            end_node = opt_to_graph[item]
                            if elem not in installed_comp:
                                comp_enabled = 1 * (yc[i] == 1)
                                installed_comp.append(elem)
                            else:
                                comp_enabled = 0
                            diam = diameter[yd[i]]
                            opt_path, opt_path_length, opt_path_cost = self.expansion(diam, comp_enabled, start_node,
                                                                                      end_node)
                            expansion_costs += opt_path_cost
                            expansion_result[(elem, item)] = [opt_path_length, diam, comp_enabled]

                    start_node = opt_to_graph[item]
                    end_node = opt_to_graph[i + 11]
                    comp_enabled = 0
                    diam = diameter[yd[i]]
                    opt_path, opt_path_length, opt_path_cost = self.expansion(diam, comp_enabled, start_node,
                                                                              end_node)
                    expansion_costs += opt_path_cost
                    expansion_result[(item, i + 11)] = [opt_path_length, diam, comp_enabled]

            else:
                for elem in start_nodes:
                    start_node = opt_to_graph[elem]
                    end_node = opt_to_graph[i + 11]
                    if elem not in installed_comp:
                        comp_enabled = 1 * (yc[i] == 1)
                        installed_comp.append(elem)
                    else:
                        comp_enabled = 0
                    diam = diameter[yd[i]]
                    opt_path, opt_path_length, opt_path_cost = self.expansion(diam, comp_enabled, start_node, end_node)
                    expansion_costs += opt_path_cost
                    expansion_result[(elem, i + 11)] = [opt_path_length, diam, comp_enabled]

        self.preparation(expansion_result)
        operation_costs = 0
        load_shedding=[918.4e3*365/2*5*yp[0], 1.5501e6*365/2*10*yp[1], 1.1921e6*365/2*5*yp[2]]
        lst=sum(load_shedding)

        for i in range(10):
            demand_factor = 1 + (i / 2 * (i % 2 == 0) + math.floor(i / 2) * (i % 2 != 0)) / 100
            demand_factor = demand_factor * (1 * (i % 2 == 0) + 1.05 * (i % 2 != 0))
            operation_costs += self.operation(demand_factor) * 365 / 2

        self.reset_assign()

        # orm = math.floor(math.log(operation_costs, 10))
        # obj = expansion_costs + 10 ** (9 - orm) * operation_costs
        obj = expansion_costs + operation_costs+lst
        self.counter += 1
        if self.counter > 100:
            self.fevalcount += 1
            self.counter = 0
            self.MATeng.quit()
            print('-------------------obj func calc count----------------------')
            print(self.fevalcount)
            self.MATeng = me.start_matlab()
        return obj


exp_test = my_minlp()


def callback_de(xk, convergence):
    if all(val >= 0 for val in exp_test.exp_cons(xk)):
        obj_list.append(exp_test.fitness(xk))
    else:
        obj_list.append(10 ** 10)

    conv_list.append(convergence)


con = lambda x: exp_test.exp_cons(x)
cons = NonlinearConstraint(con, 0, np.inf)
obj_func = lambda x: exp_test.fitness(x)
# callb = lambda x: callback_de(x)
lb = np.zeros(nT, dtype=int)
# lb = np.zeros(nT+3, dtype=int)

ub_aux = np.multiply(np.ones(nT, dtype=int), 1)
ub_aux[0] = 3
ub_aux[2] = ub_aux[6] = ub_aux[10] = 2
ub_aux[4] = 2
ub_aux[8] = 4

ub = ub_aux

# ub3 = np.multiply(np.ones(3, dtype=int), 2)
# ub = np.append(ub1, ub3)
bounds = Bounds(lb, ub)
int_cons = np.ones(nT, dtype=bool)
tic = time.time()
# xinit=[1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0]
result = differential_evolution(obj_func, bounds, disp=True, constraints=cons,callback=callback_de,
                                integrality=int_cons, maxiter=400)
print(time.time() - tic)
print(result.x)
print(result)
with open(r'obj_list.txt', 'w') as fp:
    fp.writelines(str(obj_list))
with open(r'conv_list.txt', 'w') as fp:
    fp.writelines(str(conv_list))


