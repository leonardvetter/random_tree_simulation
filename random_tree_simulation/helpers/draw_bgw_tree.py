"""
Draw a rooted plane tree from a depth-first sequence of children counts.

Input: a list of nonnegative integers c[0],...,c[n-1], where c[i] is the number
of children of node i, listed in *depth-first (preorder)* order. For any valid
rooted tree, sum(c) == n-1.

Two layouts are provided:
  - "tidy": classic layered tree (root at top or bottom)
  - "radial": root at center, branches outward with angles respecting the plane order
  - "spring": TODO

"""
from __future__ import annotations
import math
import argparse
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import networkx as nx
import scipy as sp
# ----------------------------- Tree reconstruction ----------------------------- #

def reconstruct_plane_tree(children_counts: List[int]) -> Dict[int, List[int]]:
    """
    Rebuild the rooted plane tree (children order preserved) from preorder counts.
    Returns adjacency as {parent: [child0, child1, ...]} with root=0 always present.
    """
    n = len(children_counts)
    children: Dict[int, List[int]] = {i: [] for i in range(n)}

    stack: List[Tuple[int, int]] = []  # (node_id, remaining_children_to_attach)
    for i, k in enumerate(children_counts):
        # attach to current parent (top of stack) if any
        if stack:
            parent, rem = stack[-1]
            children[parent].append(i)
            rem -= 1
            stack[-1] = (parent, rem)

        # push current node with its number of children
        stack.append((i, k))

        # pop any finished nodes
        while stack and stack[-1][1] == 0:
            stack.pop()

    if stack:
        # Should be empty if the sequence was consistent
        raise ValueError("Sequence ended prematurely (unmatched remaining children).")
    return children


# ----------------------------- Layout algorithms ----------------------------- #

def tidy_layout(children: Dict[int, List[int]],
                root: int = 0,
                ygap: float = 1.0) -> Dict[int, Tuple[float, float]]:
    """
    Reingold–Tilford-style simple tidy layout without shifts:
    leaves are placed at consecutive integer x; internal nodes at mean(child x).
    Returns positions {node: (x, y)} with y = depth (root at y=0).
    """
    x: Dict[int, float] = {}
    y: Dict[int, float] = {}
    next_x = [0.0]  # mutable closure

    def dfs(u: int, depth: int) -> None:
        y[u] = depth * ygap
        if not children[u]:  # leaf
            x[u] = next_x[0]
            next_x[0] += 1.0
        else:
            for v in children[u]:
                dfs(v, depth + 1)
            x[u] = sum(x[v] for v in children[u]) / len(children[u])

    dfs(root, 0)
    return {u: (x[u], y[u]) for u in children}


def radial_layout(children: Dict[int, List[int]],
                  root: int = 0,
                  rstep: float = 1.0,
                  angle_span: Tuple[float, float] = (0.0, 2 * math.pi)
                  ) -> Dict[int, Tuple[float, float]]:
    """
    Radial "balloon" layout:
      - each leaf occupies a small angular slot in [angle_span],
      - each internal node sits at the average angle of its descendants,
      - radius = depth * rstep.
    Preserves the plane order (left-to-right -> increasing angle).
    """
    angle: Dict[int, float] = {}
    radius: Dict[int, float] = {}

    # first get leaf order counts
    def count_leaves(u: int) -> int:
        if not children[u]:
            return 1
        return sum(count_leaves(v) for v in children[u])

    total_leaves = count_leaves(root)
    theta0, theta1 = angle_span
    slot = (theta1 - theta0) / max(1, total_leaves)
    leaf_cursor = [0]

    def dfs(u: int, depth: int) -> None:
        radius[u] = depth * rstep
        if not children[u]:
            angle[u] = theta0 + (leaf_cursor[0] + 0.5) * slot
            leaf_cursor[0] += 1
        else:
            for v in children[u]:
                dfs(v, depth + 1)
            angle[u] = sum(angle[v] for v in children[u]) / len(children[u])

    dfs(root, 0)
    pos = {u: (radius[u] * math.cos(angle[u]),
               radius[u] * math.sin(angle[u])) for u in children}
    return pos

def to_networkx(children: Dict[int, List[int]]) -> nx.Graph:
    """Convert adjacency map to an undirected NetworkX graph for layout."""
    G = nx.Graph()
    G.add_nodes_from(children.keys())
    for u, kids in children.items():
        for v in kids:
            G.add_edge(u, v)
    return G

def spring_layout_nx(children: Dict[int, List[int]],
                     root: int = 0,
                     seed: int | None = None,
                     k: float | None = None,
                     iterations: int = 1000,
                     scale: float = 1.0,
                     init: str = "radial",            # "none" | "tidy" | "radial"
                     rstep: float = 1.0,
                     fix_root: bool = False) -> Dict[int, Tuple[float, float]]:
    """
    Use networkx.spring_layout to place nodes.
      - Optionally initialize from 'tidy' or 'radial' to preserve some structure.
      - Optionally fix the root at the origin while the rest relax.
      - k controls the optimal edge length (None lets NX choose ~1/sqrt(n)).
    """
    G = to_networkx(children)

    # Optional initialization
    pos0 = None
    fixed = None
    if init.lower() == "tidy":
        pos0 = tidy_layout(children, root=root, ygap=1.0)
    elif init.lower() == "radial":
        pos0 = radial_layout(children, root=root, rstep=rstep)

    if pos0 is not None:
        # Ensure floats and (optionally) put root at center for a nice anchor
        pos0 = {u: (float(x), float(y)) for u, (x, y) in pos0.items()}
        if fix_root:
            pos0[root] = (0.0, 0.0)
            fixed = [root]
    elif fix_root:
        # If no init but root must be fixed, give a trivial init putting root at origin
        pos0 = {u: (0.0, 0.0) for u in children}
        fixed = [root]

    pos = nx.spring_layout(
        G,
        pos=pos0,           # initial guess (None -> random)
        fixed=fixed,        # nodes whose positions are fixed
        seed=seed,          # reproducible randomness
        k=k,                # optimal distance between nodes
        iterations=iterations,
        scale=scale,
        center=(0.0, 0.0),
        weight=None,        # unweighted tree edges
        threshold=1e-4,
        dim=2
    )

    # networkx returns a dict of numpy arrays; convert to plain floats
    return {u: (float(p[0]), float(p[1])) for u, p in pos.items()}


# ----------------------------- Drawing ----------------------------- #

def draw_tree(children: Dict[int, List[int]],
              pos: Dict[int, Tuple[float, float]],
              out: str = "tree.png",
              node_size: float = 6.0,
              line_w: float = 0.6,
              labels: bool = False,
              margin_scale: float = 0.08,
              invert_y: bool = False,
              figsize_scale: float = 0.28) -> None:
    """
    Render the tree with minimalist styling (white background, thin black lines, small dots).
    """
    xs = [pos[u][0] for u in pos]
    ys = [pos[u][1] for u in pos]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = xmax - xmin
    dy = ymax - ymin
    # pick a figure size that scales with the tree footprint (kept modest)
    w = max(4.0, figsize_scale * max(1.0, dx + 1))
    h = max(3.0, figsize_scale * max(1.0, dy + 1))
    fig, ax = plt.subplots(figsize=(w, h), dpi=300)

    # edges
    for u, kids in children.items():
        x0, y0 = pos[u]
        for v in kids:
            x1, y1 = pos[v]
            ax.plot([x0, x1], [y0, y1], lw=line_w, color="black", solid_capstyle="round")

    # nodes
    ax.scatter(xs, ys, s=node_size, c="black", zorder=3)

    # optional labels
    if labels:
        for u in children:
            ax.text(pos[u][0], pos[u][1], str(u), fontsize=6,
                    ha="center", va="center", color="black",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))

    # aesthetics
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    # margins
    mx = dx * margin_scale if dx > 0 else 0.5
    my = dy * margin_scale if dy > 0 else 0.5
    ax.set_xlim(xmin - mx, xmax + mx)
    ax.set_ylim(ymin - my, ymax + my)
    if invert_y:
        ax.invert_yaxis()

    plt.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)