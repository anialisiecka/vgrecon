#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void print128(unsigned __int128 n) {
    if (n == 0) {
        std::cout << 0 << "\n";
        return;
    }
    std::string s;
    while (n > 0) {
        s += '0' + (n % 10);
        n /= 10;
    }
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        std::cout << *it;
    }
    std::cout << "\n";
}

int main(int argc, char* argv[]) {

    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        std::cerr << "Cannot open the file: " << argv[1] << "\n";
        return 1;
    }

    // Run 1
    size_t max_vtx_id = 0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == 'S') {
            const char* ptr = line.data();
            const char* end = ptr + line.length();

            while (ptr < end && *ptr != '\t') ptr++;
            if (ptr >= end) continue;
            ptr++;

            size_t vtx_id = 0;
            while (ptr < end && *ptr >= '0' && *ptr <= '9') {
                vtx_id = vtx_id * 10 + (*ptr - '0');
                ptr++;
            }

            if (vtx_id > max_vtx_id) {
                max_vtx_id = vtx_id;
            }
        }
    }

    if (max_vtx_id == 0) {
        std::cout << "Homology relation size: 0\n";
        return 0;
    }

    std::vector<size_t> vertex_lengths(max_vtx_id + 1, 0);
    std::vector<size_t> vertex_counts(max_vtx_id + 1, 0);

    file.clear();
    file.seekg(0, std::ios::beg);

    // Run 2
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == 'S') {
            const char* ptr = line.data();
            const char* end = ptr + line.length();

            while (ptr < end && *ptr != '\t') ptr++;
            if (ptr >= end) continue;
            ptr++;

            size_t vtx_id = 0;
            while (ptr < end && *ptr >= '0' && *ptr <= '9') {
                vtx_id = vtx_id * 10 + (*ptr - '0');
                ptr++;
            }

            while (ptr < end && *ptr != '\t') ptr++;
            if (ptr >= end) continue;
            ptr++;

            const char* seq_start = ptr;
            while (ptr < end && *ptr != '\t' && *ptr != '\r' && *ptr != '\n') {
                ptr++;
            }
            vertex_lengths[vtx_id] = ptr - seq_start;
        }
        else if (line[0] == 'P') {
            size_t t1 = line.find('\t');
            size_t t2 = line.find('\t', t1 + 1);
            size_t t3 = line.find('\t', t2 + 1);
            if (t2 == std::string::npos) continue;
            if (t3 == std::string::npos) t3 = line.length();

            const char* path_ptr = line.data() + t2 + 1;
            const char* path_end = line.data() + t3;

            while (path_ptr < path_end) {
                const char* comma = path_ptr;
                while (comma < path_end && *comma != ',') {
                    comma++;
                }

                if (comma - path_ptr < 2) break;

                size_t vtx_id = 0;
                for (const char* p = path_ptr; p < comma - 1; ++p) {
                    vtx_id = vtx_id * 10 + (*p - '0');
                }

                vertex_counts[vtx_id]++;

                path_ptr = comma + 1;
            }
        }
    }
    file.close();

    // final homology relation calculation
    unsigned __int128 total_sum = 0;

    for (size_t vtx_id = 1; vtx_id <= max_vtx_id; ++vtx_id) {
        size_t count = vertex_counts[vtx_id];
        if (count < 2) continue;

        size_t len = vertex_lengths[vtx_id];

        unsigned __int128 combinations = (static_cast<unsigned __int128>(count) * (count - 1)) / 2;
        total_sum += static_cast<unsigned __int128>(len) * combinations;
    }

    std::cout << "Homology relation size: ";
    print128(total_sum);

    return 0;
}
