//
// Created by fritsche on 05/04/2022.
//

#pragma once

#include <cstdint>
#include <vector>
#include <stack>

/**
 * ttree is a custom n-ary tree library
 */
namespace ttree {
    // Helper function that allows for static_assert(false, ...)
    // "false" must be always_false<t> and thus depend on a type to not be
    // immediately caught out by the compiler
    template <class... T>
    constexpr bool always_false = false;

    auto comp_int = [](int a, int b) {
        return a < b;
    };

    template <typename T>
    struct InsertPolicy {
        // returns pair of Inorder id and level
        std::pair<uint64_t, int64_t> GetTreePosition (T& element) const {
            static_assert(always_false<T>, "Please provide your own Insert Policy for the provided template type");
            return { 0, -1 };
        }
    };

    template<>
    struct InsertPolicy <int> {
        // returns pair of Inorder id and level
        std::pair<uint64_t, int64_t> GetTreePosition (int& element) const {
            return { element, -1 };
        }
    };

    static constexpr int AnswerToEverything() {
        return 42;
    }

    using IdxType = uint32_t;
    template<typename T, class Compare, class InsertPolicy>
    struct node : private InsertPolicy {
        T element;
        IdxType node_index = 0;
        IdxType parent_index = 0;

        std::vector<IdxType> children;

        node(T& element, IdxType node_index) : element(element), node_index(node_index) {};
        node(T& element, IdxType node_index, IdxType parent_index) : element(element), node_index(node_index), parent_index(parent_index) {};

        const std::string ToString() const {
            std::string result = "(element, index, parent_index) ";
            result += std::to_string(element);
            result += ", ";
            result += std::to_string(node_index);
            result += ", ";
            result += std::to_string(parent_index);
            return result;
        }

        const bool IsRoot() const {
            return node_index == parent_index;
        }

        inline const std::vector<IdxType>& Children() const {
            return children;
        }

        inline void AddChild(const IdxType& child_index) {
            children.emplace_back(child_index);
            std::sort(children.begin(), children.end(), Compare());
        }
    };



    template<typename T, class Compare, class InsertPolicy>
    class TinyTree : private InsertPolicy {
    public:
        using Node = node<T, Compare, InsertPolicy>;

    private:
        IdxType m_root_index = 0;
        std::vector<Node> m_nodes;

        void PrintTreeWorker(IdxType node_index, int32_t depth) {
            assert(node_index < m_nodes.size() && "node_index is invalid.");
            std::cout << m_nodes[node_index].ToString() << std::endl;
            depth++;

            for (auto& child_idx : m_nodes.Children()) {
                PrintTreeWorker(child_idx, depth);
            }
        }

        bool IsValidIdx(IdxType idx) {
            return idx >= 0 && idx < m_nodes.size();
        }

    public:
        using insert_policy = InsertPolicy;

        TinyTree() {}

        Node& Root() {
            return m_nodes[m_root_index];
        }

        inline constexpr size_t Size() const {
            return m_nodes.size();
        }

        inline uint64_t GetInorderId(T& value) {
            return this->GetTreePosition(value).first;
        }

        template<typename Functor>
        void InorderTraversal(Functor functor, Node& start_node) {

            // Explicitly signed
            std::cout << "start with index " << start_node.node_index << std::endl;

            std::stack<IdxType> stack;
            stack.push(start_node.node_index);

            while (!stack.empty()) {
                auto& current_node = m_nodes[stack.pop()];
                functor(current_node);

                auto begin = current_node.Children().end();
                while (begin-- != current_node.Children().begin()) {
                    stack.push(*begin);
                }
            }
        }

        template<typename Functor>
        void InorderTraversal(Functor functor) {
            InorderTraversal(functor, m_nodes[m_root_index]);
        }

        void PrintTree() {
            if (!m_nodes.empty()) {
                PrintTreeWorker(m_root_index,0);
            }
        }

        // Return pair of
        const std::pair<int64_t, IdxType> FindIdx(const T& value) const {
            if (m_nodes.empty())
                return { -1, m_root_index };

            const auto [ target_inorder_id, target_level ] = this->GetTreePosition(value);

            auto& current_node = Root();

            uint64_t current_inorder_id = 0;
            int64_t current_level = 0;

            uint64_t child_inorder_id = 0;
            int64_t child_level = 0;

            std::tie(current_inorder_id, current_level) = this->GetTreePosition(current_node.element);

            while (true) {
                // Target level is larger or equal to the current node level
                // Target inorder_id is larger than the current node
                // This means the target node could be
                assert(std::is_sorted(current_node.Children().begin(), current_node.Children().end(), [this](IdxType &a, IdxType &b) {
                    return this->GetTreePosition(m_nodes[a].element).first < this->GetTreePosition(m_nodes[b].element).first;
                }));


                auto last_child_idx = current_node.node_idx;
                for (auto& child_idx : current_node.Children()) {
                    std::tie(child_inorder_id, child_level) = this->GetTreePosition(m_nodes[child_idx]);

                    if (target_inorder_id > child_inorder_id) {
                        // Either target is in subtree of last_child_idx
                        // Or target should replace last_child_idx

                        std::tie(child_inorder_id, child_level) = this->GetTreePosition(m_nodes[last_child_idx]);
                        if (child_level > target_level) {
                            // Element should be located as child of current_node but replace the child last_child_idx
                            return { last_child_idx, current_node.node_idx };
                        } else {
                            // search on in child
                            current_node = m_nodes[last_child_idx];
                            break;
                        }
                    }
                    last_child_idx = child_idx;
                }

                std::tie(current_inorder_id, current_level) = this->GetTreePosition(current_node.element);
            }

            return { 0, 0 };
        }

        bool Contains(T& value) {
            const auto [ target_inorder_id, target_level ] = this->GetTreePosition(value);
            auto [ target_index, parent_index ] = FindIdx(target_inorder_id);

//            std::cout << "Contains " << value << "? (target_index, parent_index): " << target_index << ", " << parent_index << std::endl;

            return target_index != -1;
        }


        const IdxType Insert(T value) {
            if (m_nodes.empty()) {
                m_nodes.emplace_back( Node{ value, 0, 0 } );
                return m_root_index;
            }

            const auto [ target_inorder_id, target_level ] = this->GetTreePosition(value);

            if (target_level == -1) {
                auto [ target_index, parent_index ] = FindIdx(target_inorder_id);
//                std::cout << "Ins: (target_index, parent_index): " << target_index << ", " << parent_index << std::endl;

                if (target_index != -1) {
                    return target_index;
                }

                IdxType insert_index = m_nodes.size();
                m_nodes.emplace_back( Node{ value, insert_index, parent_index } );
                auto& parent_node = m_nodes[parent_index];


                const auto [ parent_inorder_id, parent_level ] = this->GetTreePosition(parent_node.element);

                assert(target_inorder_id != parent_inorder_id && "Trying to insert value that is already in tree.");

                if (target_inorder_id < parent_inorder_id)
                    parent_node.SetLeftChild(insert_index);
                if (target_inorder_id > parent_inorder_id)
                    parent_node.SetRightChild(insert_index);


                return m_nodes.size() - 1;

            } else {
                // To be implemented
            }
            return 0;
        }
    };
}