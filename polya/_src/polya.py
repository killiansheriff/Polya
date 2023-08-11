import itertools

import igraph
import numpy as np
import sympy as sp
from sympy.combinatorics.perm_groups import PermutationGroup


class Polya:
    """Class to get Polya's pattern inventory.

    ```python
    polya = Polya(graph_name="fcc", ntypes=3)
    poly, inms = polya.get_gt()
    ```
    """

    def __init__(self, graph_name, ntypes):
        """Class to get Polya's pattern inventory.

        ```python
        polya = Polya(graph_name="fcc", ntypes=3)
        poly, inms = polya.get_gt()
        ```
        """
        self.graph_generators = {
            "fcc": self.get_fcc_1nn_graph,
            "bcc": self.get_bcc_1nn_graph,
            "hcp": self.get_hcp_1nn_graph,
            "fcc_1nn2nn3nn": self.get_fcc_1nn_2nn_3nn_graph,
        }
        self.graph_name = graph_name
        self.ntypes = ntypes
        self.graph_generator = self.graph_generators[self.graph_name]

    def get_edges(self, vertexpositions, nn_dst, atol=0.1):
        # Subtract each point from all the other points
        diff = vertexpositions[:, np.newaxis, :] - vertexpositions[np.newaxis, :, :]
        # Compute the norm of the differences along the last axis
        distances = np.linalg.norm(diff, axis=-1)
        # Add edges to the graph
        mask = np.isclose(distances, nn_dst, atol=atol)
        edge_index_0, edge_index_1 = np.where(mask)
        edges = [(edge_index_0[i], edge_index_1[i]) for i in range(len(edge_index_0))]
        edges = np.unique(np.sort(edges), axis=0)
        return edges, distances

    def generate_permutations(self, values, coords):
        signs = list(itertools.product([1, -1], repeat=len(values)))
        permutations = []
        for sign in signs:
            value = tuple(val * s for val, s in zip(values, sign))
            permutations.extend(list((itertools.permutations(value, 3))))
        coords = np.vstack((coords, np.array(list(set(permutations)))))

        return coords

    def get_fcc_1nn_2nn_3nn_graph(self):
        # Central atom (Origin)
        coords = np.zeros((1, 3))
        # 1NN atoms (sqrt(2)/2 away)
        values = (0, 1 / 2, 1 / 2)
        coords = self.generate_permutations(values, coords)

        # 2NN atoms (1.0 away)
        values = (0, 0, 1)
        coords = self.generate_permutations(values, coords)

        # 3NN atoms (sqrt(3/2) away)
        values = (1, 1 / 2, 1 / 2)
        coords = self.generate_permutations(values, coords)

        # remove central atom
        vertexpositions = coords[1:]

        g = igraph.Graph()
        g.add_vertices(len(vertexpositions))
        edges, distances = self.get_edges(vertexpositions, np.sqrt(2) / 2)
        g.add_edges(edges)
        g.vs["pos"] = vertexpositions

        # Quick viz plotting
        # import matplotlib.pyplot as plt
        # fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        # h = g
        # # igraph draw
        # ax1.set_title("Plot with igraph plot")
        # layout = h.layout_kamada_kawai()
        # igraph.plot(h, layout=layout, target=ax1)
        # plt.axis("off")
        # plt.savefig("1nn2nn3nn.png")

        return g

    def get_fcc_1nn_graph(self):
        # Define the vertex positions
        C0 = 0.816496580927726032732428024902
        vertexpositions = np.array(
            [
                # (0, 0, 0),
                (C0, 0.0, C0),
                (C0, 0.0, -C0),
                (-C0, 0.0, C0),
                (-C0, 0.0, -C0),
                (C0, C0, 0.0),
                (C0, -C0, 0.0),
                (-C0, C0, 0.0),
                (-C0, -C0, 0.0),
                (0.0, C0, C0),
                (0.0, C0, -C0),
                (0.0, -C0, C0),
                (0.0, -C0, -C0),
            ]
        )

        g = igraph.Graph()
        g.add_vertices(len(vertexpositions))
        edges, distances = self.get_edges(vertexpositions, 1.154)
        g.add_edges(edges)
        g.vs["pos"] = vertexpositions

        return g

    def get_fcc_1nn2nn_graph(self):
        # Define the vertex positions
        C0 = 0.816496580927726032732428024902
        C1 = 1.035276180410083049395595350496

        vertexposition = [
            (0.0, 0.0, C1),
            (0.0, 0.0, -C1),
            (C1, 0.0, 0.0),
            (-C1, 0.0, 0.0),
            (0.0, C1, 0.0),
            (0.0, -C1, 0.0),
            (C0, 0.0, C0),
            (C0, 0.0, -C0),
            (-C0, 0.0, C0),
            (-C0, 0.0, -C0),
            (C0, C0, 0.0),
            (C0, -C0, 0.0),
            (-C0, C0, 0.0),
            (-C0, -C0, 0.0),
            (0.0, C0, C0),
            (0.0, C0, -C0),
            (0.0, -C0, C0),
            (0.0, -C0, -C0),
        ]

        # Create an empty igraph graph object
        g = igraph.Graph()

        # Add vertices to the graph
        g.add_vertices(len(vertexposition))

        # Add edges to the graph
        edges = (
            np.array(
                [
                    (15, 9),
                    (9, 17),
                    (17, 7),
                    (7, 15),
                    (10, 18),
                    (18, 8),
                    (8, 16),
                    (16, 10),
                    (8, 12),
                    (8, 11),
                    (11, 7),
                    (7, 12),
                    (9, 14),
                    (14, 10),
                    (10, 13),
                    (13, 9),
                    (15, 13),
                    (13, 16),
                    (16, 11),
                    (11, 15),
                    (14, 17),
                    (17, 12),
                    (12, 18),
                    (18, 14),
                    (1, 15),
                    (1, 9),
                    (1, 17),
                    (1, 7),
                    (5, 11),
                    (5, 15),
                    (5, 13),
                    (5, 16),
                    (2, 8),
                    (2, 16),
                    (2, 18),
                    (2, 10),
                    (3, 7),
                    (3, 11),
                    (3, 8),
                    (3, 12),
                    (6, 17),
                    (6, 14),
                    (6, 18),
                    (6, 12),
                    (4, 9),
                    (4, 13),
                    (4, 10),
                    (4, 14),
                ]
            )
            - 1
        )

        g.add_edges(edges)

        # Set the vertex positions
        g.vs["pos"] = vertexposition

        return g

    def get_bcc_1nn_graph(self):
        vertexpositions = np.array(
            [
                # [0, 0, 0],
                [1, 1, 1],
                [1, 1, -1],
                [1, -1, 1],
                [1, -1, -1],
                [-1, 1, 1],
                [-1, 1, -1],
                [-1, -1, 1],
                [-1, -1, -1],
            ]
        )
        # Create an empty igraph graph object
        g = igraph.Graph(directed=False)

        g.add_vertices(len(vertexpositions))

        edges_1nn, distances = self.get_edges(vertexpositions, 1.73)
        edges_2nn, _ = self.get_edges(vertexpositions, 2)

        if len(edges_1nn) > 0:
            edges = np.concatenate([edges_1nn, edges_2nn])
        else:
            edges = edges_2nn

        g.add_edges(edges)

        return g

    def get_hcp_1nn_graph(self):
        # Define positions of nearest neighbors relative to central atom at (0, 0, 0)

        # Define the lattice parameters
        a = 1
        c = np.sqrt(8 / 3)

        # Define the positions of the first nearest neighbors
        # vertexpositions = np.array(
        #     [
        #         (0, 0, 0),
        #         (1, 0, 0),
        #         (-1, 0, 0),

        #         (1 / 2, np.sqrt(3) * a / 2, 0),
        #         (-a / 2, -np.sqrt(3) * a / 2, 0),

        #         (a / 2, 0, c / 2),
        #         (-a / 2, 0, c / 2),
        #         (0, np.sqrt(3) * a / 2, c / 2),
        #         (0, -np.sqrt(3) * a / 2, c / 2),

        #         (a, np.sqrt(3) * a / 3, c / 2),
        #         (-a, -np.sqrt(3) * a / 3, c / 2),
        #         (a, -np.sqrt(3) * a / 3, c / 2),
        #         (-a, np.sqrt(3) * a / 3, c / 2),
        #     ]
        # )

        vertexpositions = np.array(
            [
                # [0, 0, 0],
                [1, 0, 0],
                [-1, 0, 0],
                [1 / 2, np.sqrt(3) / 2, 0],
                [1 / 2, -np.sqrt(3) / 2, 0],
                [-1 / 2, np.sqrt(3) / 2, 0],
                [-1 / 2, -np.sqrt(3) / 2, 0],
                [0, np.sqrt(3) / 3, 1 / 2],
                [0, np.sqrt(3) / 3, -1 / 2],
                [1 / 2, -np.sqrt(3) / 6, 1 / 2],
                [1 / 2, -np.sqrt(3) / 6, -1 / 2],
                [-1 / 2, -np.sqrt(3) / 6, 1 / 2],
                [-1 / 2, -np.sqrt(3) / 6, -1 / 2],
            ]
        )
        vertexpositions *= np.array([a, a, c])

        # Create an empty igraph graph object
        g = igraph.Graph(directed=False)

        # Add vertices to the graph
        g.add_vertices(len(vertexpositions))
        edges, distances = self.get_edges(vertexpositions, 1, atol=0.1)
        g.add_edges(edges)

        # Set the vertex positions
        g.vs["pos"] = vertexpositions

        return g

    def get_cycle_index(self, permgroup):
        cycle_types = [p.cycle_structure for p in permgroup]
        monomials = [
            np.prod([sp.symbols(f"s_{ctype}") ** cycle[ctype] for ctype in cycle.keys()])
            for cycle in cycle_types
        ]
        nnodes = np.sum([key * value for key, value in cycle_types[0].items()])
        group_size = len(permgroup) + 1  # add identity
        cycle_index = np.sum(monomials) + sp.symbols(f"s_1") ** nnodes
        return cycle_index / group_size  # need divided size of group

    def get_gt(self):
        self.g = self.graph_generator()

        nnodes = self.g.vcount()

        # Compute the automorphism group
        permgroup = np.array(self.g.get_automorphisms_vf2())

        # Get the permutation representation of the group
        permgroup = PermutationGroup(permgroup)
        cycle_index = self.get_cycle_index(permgroup)

        # define symbolic variables for d1 to d10
        types = sp.symbols(f"t1:{self.ntypes+1}")

        # replace s_i with the sum of the powers of the d variables and factorize
        dpoly = sp.factor(
            cycle_index.subs(
                [
                    (sp.symbols(f"s_{i}"), np.sum([types[j] ** i for j in range(self.ntypes)]))
                    for i in range(1, nnodes + 1)
                ]
            )
        )

        # replace s_i with the sum of the powers of 1 for each variable and factorize
        onepoly = sp.factor(
            cycle_index.subs(
                [
                    (sp.symbols(f"s_{i}"), sum([1**i for _ in range(self.ntypes)]))
                    for i in range(1, nnodes + 1)
                ]
            )
        )
        return dpoly, onepoly
