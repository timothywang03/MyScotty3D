#include "test.h"
#include "geometry/halfedge.h"

static void expect_flip(Halfedge_Mesh &mesh, Halfedge_Mesh::EdgeRef edge, Halfedge_Mesh const &after) {
	if (auto ret = mesh.flip_edge(edge)) {
		if (auto msg = mesh.validate()) {
			throw Test::error("Invalid mesh: " + msg.value().second);
		}
		// check if returned edge is the same edge
		if (ret != edge) {
			throw Test::error("Did not return the same edge!");
		}
		// check mesh shape:
		if (auto difference = Test::differs(mesh, after, Test::CheckAllBits)) {
			throw Test::error("Resulting mesh did not match expected: " + *difference);
		}
	} else {
		throw Test::error("flip_edge rejected operation!");
	}
}

/*
BASIC CASE

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Flip Edge on Edge: 1-4

After mesh:
0--1\
|\   \
| \---2
|    /
3--4/
*/
Test test_a2_l1_flip_edge_basic_simple("a2.l1.flip_edge.basic.simple", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 2}, 
		{0, 2, 1}
	});

	expect_flip(mesh, edge, after);
});

/*
EDGE CASE

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Flip Edge on Edge: 3-4

After mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/
*/
Test test_a2_l1_flip_edge_edge_boundary("a2.l1.flip_edge.edge.boundary", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

	if (mesh.flip_edge(edge)) {
		throw Test::error("flip_edge should not work at the boundary.");
	}
});

/*
MULTI-FLIP CASE

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Flip Edge on Edge: 1-4->0-2->3-1

After1 mesh:
0--1\
|\   \
| \---2
|    /
3--4/

After2 mesh:
0--1\
|  | \
|  /  2
|/   /
3--4/

After3 mesh:
0--1\
||   \
| |   2
|  | /
3--4/

After4 mesh:
0--1\
|    \
| ----2
|/   /
3--4/

After5 mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/
*/
Test test_a2_l1_flip_edge_basic_multi("a2.l1.flip_edge.basic.multi", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {{0, 3, 4, 1}, {1, 4, 2}});
  Halfedge_Mesh::EdgeRef edge1 = mesh.halfedges.begin()->next->next->edge;

  Halfedge_Mesh after1 = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {{0, 3, 4, 2}, {0, 2, 1}});

  Halfedge_Mesh::EdgeRef edge2 =
      after1.halfedges.begin()->next->next->next->edge;

  Halfedge_Mesh after2 = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {{3, 4, 2, 1}, {1, 0, 3}});

  Halfedge_Mesh::EdgeRef edge3 =
      after2.halfedges.begin()->next->next->next->edge;

  Halfedge_Mesh after3 = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {
          {0, 4, 2, 1},
          {0, 3, 4},
      });

  Halfedge_Mesh::EdgeRef edge4 = after3.halfedges.begin()->edge;

  Halfedge_Mesh after4 = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {
          {2, 1, 0, 3},
          {2, 3, 4},
      });

  Halfedge_Mesh::EdgeRef edge5 =
      after4.halfedges.begin()->next->next->next->edge;

  Halfedge_Mesh last = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(2.2f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {{0, 3, 4, 1}, {1, 4, 2}});

  expect_flip(mesh, edge1, after1);
  printf("\nPassed first flip!\n");
  expect_flip(after1, edge2, after2);
  printf("Passed second flip!\n");
  expect_flip(after2, edge3, after3);
  printf("Passed third flip!!\n");
  expect_flip(after3, edge4, after4);
  printf("Passed fourth flip!!\n");
  expect_flip(after4, edge5, last);
  printf("Passed final flip!!!\n");
});

/*
Edge CASE

Initial mesh:
0--1--2
|  |  |
|  |  |
|  6  |
|  |  |
|  |  |
3--4--5

Flip Edge on Edge: 1-5

After mesh:
0--1--2
|  |  |
|  |  |
|  6  |
|  |  |
|  |  |
3--4--5
*/
Test test_a2_l1_flip_edge_edge_case("a2.l1.flip_edge.case", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f), Vec3(4.0f, 0.0f, 0.0f),
       Vec3(0.0f, 4.0f, 0.0f), Vec3(2.0f, 4.0f, 0.0f), Vec3(4.0f, 4.0f, 0.0f),
       Vec3(2.0f, 2.0f, 0.0f)},
      {{0, 1, 6, 4, 3}, {1, 2, 5, 4, 6}});
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

  if (mesh.flip_edge(edge)) {
    throw Test::error(
        "flip_edge should not work creates wrong mesh - square case.");
  }
});

/*
Edge CASE

Initial mesh:
0--1--2
|  \  |
|    7|
|  /  |
| 6   |
|  \  |
3--4--5

Flip Edge on Edge: 6-7

After mesh:
0--1--2
|  \  |
|    7|
|  /  |
| 6   |
|  \  |
3--4--5
*/
Test test_a2_l1_flip_edge_zigzag_edge_case("a2.l1.flip_edge.case.zigzag", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {
          Vec3(0.0f, 0.0f, 0.0f),
          Vec3(2.0f, 0.0f, 0.0f),
          Vec3(4.0f, 0.0f, 0.0f),
          Vec3(0.0f, 4.0f, 0.0f),
          Vec3(2.0f, 4.0f, 0.0f),
          Vec3(4.0f, 4.0f, 0.0f),
          Vec3(1.0f, 1.0f, 0.0f),
          Vec3(3.0f, 3.0f, 0.0f),
      },
      {{0, 1, 7, 6, 4, 3}, {1, 2, 5, 4, 6, 7}});
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

  if (mesh.flip_edge(edge)) {
    throw Test::error(
        "flip_edge should not work creates wrong mesh - zigzag case.");
  }
});

/*
EDGE CASE

Initial mesh:
0---1
|  /|
| 2 |
|/  |
3---4

Flip Edge on Edge: 2-1

After mesh:
0---1
|  /|
| 2 |
|/  |
3---4
*/
Test test_a2_l1_flip_edge_edge_reject("a2.l1.flip_edge.edge.reject", []() {
  Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces(
      {Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f),
       Vec3(-1.3f, -0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)},
      {{2, 1, 0, 3}, {3, 4, 1, 2}});
  Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->edge;

  if (mesh.flip_edge(edge)) {
    if (auto msg = mesh.validate()) {
      throw Test::error("Invalid mesh: " + msg.value().second);
    }
    throw Test::error("flip_edge should not work as it rejects.");
  }
});

/*
Initial
1
|\
| \
4\ \
| \ \
|  \ \
0--3--2
*/

// Test test_a2_l2_flip_edge_triangle("a2.l1.flip_edge.edge.triangle", []() {
//     Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces(
//         {Vec3(0.0f, 0.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
//         Vec3(0.5f, 0.5f, 0.0f), Vec3()},
//         {{1, 0, 3}, {1, 3, 2}});

//     Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

    
// })