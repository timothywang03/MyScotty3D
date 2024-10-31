
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>
#include <functional>	
#include <iostream>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

	unsigned long n = primitives.size();
	unsigned long BUCKET_SIZE = 40;

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
	std::function<size_t(std::vector<Primitive>&, unsigned long, unsigned long)> build_bvh;
    build_bvh = [&](std::vector<Primitive>& primitives, unsigned long start, unsigned long end) -> size_t {
        float best_cost = std::numeric_limits<float>::infinity();
        size_t best_axis = 0;
        size_t best_partition = (start + end) / 2;

		BBox root_box;
		for (unsigned long i = start; i < end; i++) {
			root_box.enclose(primitives[i].bbox());
		}

        if (end - start <= max_leaf_size) {
            return BVH::new_node(root_box, start, end - start, 0, 0);
        }

        for (int axis = 0; axis < 3; axis++) {
            // sort the primitives along the i-th axis
            std::sort(primitives.begin() + start, primitives.begin() + end,
                        [axis](const Primitive& a, const Primitive& b) {
                        return a.bbox().center()[axis] < b.bbox().center()[axis];
                        });

			// std::cout << "partitioning" << std::endl;
            // create a bounding box for each primitive
            for (unsigned long i = start + 1; i < end - 1; i += BUCKET_SIZE) {
                BBox left, right;

                for (unsigned long j = start; j < i; j++) {
                    left.enclose(primitives[j].bbox());
                }
                for (unsigned long j = i; j < end; j++) {
                    right.enclose(primitives[j].bbox());
                }
                float cost = left.surface_area() * (i - start) + right.surface_area() * (end - i);
                if (cost < best_cost) {
                    best_cost = cost;
                    best_axis = axis;
                    best_partition = i;
                }
            }
        }
        // sort the primitives along the best axis
        std::sort(primitives.begin() + start, primitives.begin() + end,
                    [best_axis](const Primitive& a, const Primitive& b) {
                    return a.bbox().center()[best_axis] < b.bbox().center()[best_axis];
                    });
        size_t left_node = build_bvh(primitives, start, best_partition);
        size_t right_node = build_bvh(primitives, best_partition, end);

		return BVH::new_node(root_box, start, end - start, left_node, right_node);
    };

    build_bvh(primitives, 0, n);
	root_idx = nodes.size() - 1;
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	if (primitives.empty()) return Trace();

	// get the root bvh node
	Node root = nodes[root_idx];

	// keep track globally of the best hit
	float ray_traverse_ub = ray.dist_bounds.y;
	Trace best = Trace();
	// bool updated = false;
	float best_distance = std::numeric_limits<float>::infinity();

	// define a traversal that will return the closest intersection between ray and bounding box
	std::function<void(const Ray&, Node)> traverse;
	traverse = [&](const Ray& ray, Node node) -> void {
		if (node.is_leaf()) {
			// base case: if it is a leaf node, iterate through all the primitives in the node
			for (size_t i = node.start; i < node.start + node.size; i++) {
				Trace hit = primitives[i].hit(ray);
				best = Trace::min(best, hit);
			}
			if (best.hit) {
				best_distance = best.distance;
				// updated = true;
			}
			return;	// regardless of whether there is a hit, return the closest hit

		} else {
			Vec2 hit_times_close = Vec2{0.0, ray_traverse_ub};
			Vec2 hit_times_far = Vec2{0.0, ray_traverse_ub};

			Node close_node = nodes[node.l];
			Node far_node = nodes[node.r];
			bool hit_close = close_node.bbox.hit(ray, hit_times_close);
			bool hit_far = far_node.bbox.hit(ray, hit_times_far);

			if (hit_close && hit_far) {
				// if both nodes are hit, traverse the closest node first
				traverse(ray, close_node);
				if (best_distance < hit_times_far.x) {
					return;
				}
				traverse(ray, far_node);
				return;
			} else if (hit_close) {
				traverse(ray, close_node);
			} else if (hit_far) {
				traverse(ray, far_node);
			}
			return;
		}
	};

	traverse(ray, root);
	return best;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
