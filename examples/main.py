from polya import Polya
polya = Polya(graph_name="fcc", ntypes=3)
p_g, nms = polya.get_gt()
print(p_g)
