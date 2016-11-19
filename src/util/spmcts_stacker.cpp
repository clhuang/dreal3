/*********************************************************************
Author: Sicun Gao <sicung@mit.edu>

dReal -- Copyright (C) 2013 - 2016, the dReal Team

dReal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

dReal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with dReal. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#include "util/stacker.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "util/logging.h"

namespace dreal {

using std::vector;
using std::shared_ptr;
using std::numeric_limits;

spmcts_stacker::spmcts_node::spmcts_node(double const c, double const d, box const & b) : b(b), c(c), d(d) {}

spmcts_stacker::spmcts_node::spmcts_node(spmcts_node *parent, box const & b) : parent(parent), b(b), c(parent->c), d(parent->d) {}

void spmcts_stacker::spmcts_node::update(double x) {
    double delta = x - mean;
    visitcount++;
    mean += delta / visitcount;
    m2 += (x-mean)*delta;
}

double spmcts_stacker::spmcts_node::get_score() {
    double score = mean;
    if (parent)
        score += c * sqrt(2 * log(parent->visitcount) / visitcount);

    if (d >= 0)
        score += sqrt((m2+d)/visitcount);
    return score;
}

/**
 * Uses a multi-armed bandit to select the node to expand,
 * the score of each node is the reciprocal error plus
 */
spmcts_stacker::spmcts_stacker(vector<box> & boxes, scoped_vec<shared_ptr<constraint>> const & ctrs, double const prec, double const c, double const d) : stacker::stacker(boxes, ctrs, prec) {
    assert(boxes.size() == 1);
    //spNode node(new spmcts_node(c, d, boxes[0]));
    //root = std::move(node);
    root.reset(new spmcts_node(c, d, boxes[0]));
}

bool spmcts_stacker::playout() {
    vector<spmcts_node*> stack;

    if (last && !last->left) { // no child of this node
        last = NULL;
    }

    spmcts_node* curr = root.get();
    for (;;) {
        stack.push_back(curr);
        if (curr->left && curr->right) {
            if (curr->left->get_score() > curr->right->get_score()) {
                curr = curr->left.get();
            } else {
                curr = curr->right.get();
            }
        } else if (curr->left) {
            curr = curr->left.get();
        } else if (curr->right) {
            curr = curr->right.get();
        } else {
            break; //unexpanded leaf
        }
    }
    last = curr;
    curr->b.sample_point();
    double err = sample_err(curr->b);
    double score = 1.0/err;
    DREAL_LOG_INFO << score;
    for (auto node : stack) {
        node->update(score);
    }
    m_best_score = last->get_score();
    return err < m_prec;
}

box spmcts_stacker::pop_best() {
    return last->b;
}

double spmcts_stacker::sample_err(box const & b) {
    double total_err = 0;
    box sample = b.sample_point();
    
    for (auto ctr : m_ctrs) {
        assert(ctr->get_type() == constraint_type::Nonlinear);
        double err = ctr->eval_error(sample);
        DREAL_LOG_INFO << "playout current err: " << err << " obtained on ctr " << *ctr;
        total_err += err;
    }
    return total_err;
}

/**
 * Call this after playout to add children to the best node popped.
 */
void spmcts_stacker::push(box const & b) {
    assert(last != NULL);
    double score = 1.0/sample_err(b);
    spNode node(new spmcts_node(last, b));
    node->update(score);
    if (!last->left) {
        last->left = std::move(node);
    } else if (!last->right) {
        last->right = std::move(node);
    } else {
        assert(false); // something went wrong
    }
}

}  // namespace dreal;
