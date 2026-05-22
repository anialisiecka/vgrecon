#include <queue>
#include <iostream>
#include <cstdint>
#include <numeric>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <string>
#include <string_view>

struct FastaData {
    std::vector<size_t> startingPos;
    std::vector<std::string> idToSeq;
    std::unordered_map<std::string, int> seqToId;

    void init(const std::string& filename) {

        std::ios::sync_with_stdio(false);
        std::ifstream file(filename, std::ios::binary);
        if (!file) return;

        std::string line, currentIdStr;

        int idCounter = 0;
        size_t currentPos = 1;
        size_t prevPos = 1;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!currentIdStr.empty()) {
                    startingPos.push_back(prevPos);
                    prevPos = currentPos;

                    seqToId[currentIdStr] = idCounter++;
                    idToSeq.push_back(std::move(currentIdStr));
                }

                size_t firstSpace = line.find_first_of(" \t\r");
                if (firstSpace != std::string::npos) {
                    currentIdStr = line.substr(1, firstSpace - 1);
                } else {
                    currentIdStr = line.substr(1);
                }
            } else {
                if (line.back() == '\r') line.pop_back();
                currentPos += line.length();
            }
        }

        if (!currentIdStr.empty()) {
            startingPos.push_back(prevPos);
            startingPos.push_back(currentPos);
            seqToId[currentIdStr] = idCounter;
            idToSeq.push_back(std::move(currentIdStr));
        }
    }
};

struct GraphA {
    std::vector<size_t> starts;
    std::vector<long> orientedIds;

    void init(const std::string& gfa, FastaData& data) {
        std::ifstream file(gfa);
        if (!file.is_open()) return;

        std::string line;
        std::vector<size_t> vtx_lengths;
        vtx_lengths.push_back(0);

        struct PathInfo { std::string path; int seqId; };
        std::vector<PathInfo> paths;

        std::vector<size_t> counts(data.idToSeq.size(), 0);

        // Step 1
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == 'S') {
                size_t t1 = line.find('\t');
                size_t t2 = line.find('\t', t1 + 1);
                size_t t3 = line.find('\t', t2 + 1);

                // sequence length
                size_t seq_end = (t3 == std::string::npos) ? line.length() : t3;
                vtx_lengths.push_back(seq_end - (t2 + 1));

            } else if (line[0] == 'P') {
                size_t t1 = line.find('\t');
                size_t t2 = line.find('\t', t1 + 1);
                size_t t3 = line.find('\t', t2 + 1);
                if (t3 == std::string::npos) t3 = line.length();

                std::string genome = line.substr(t1 + 1, t2 - t1 - 1);
                std::string path_str = line.substr(t2 + 1, t3 - t2 - 1);

                int seqId = data.seqToId[genome];

                size_t elements = 1;
                for (char c : path_str) if (c == ',') elements++;

                counts[seqId] = elements;
                paths.push_back({std::move(path_str), seqId});
            }
        }

        // Step 2
        size_t total_elements = 0;
        std::vector<size_t> offsets(data.idToSeq.size() + 1, 0);
        for (size_t i = 0; i < data.idToSeq.size(); ++i) {
             size_t n = (counts[i] > 0) ? counts[i] : 1;
            offsets[i] = total_elements;
            total_elements += n;
        }
        offsets[data.idToSeq.size()] = total_elements;

        starts.resize(total_elements);
        orientedIds.resize(total_elements);

        // Step 3
        for (auto& p : paths) {
            size_t out_idx = offsets[p.seqId];
            size_t shift = data.startingPos[p.seqId];
            size_t current_pos = 0;

            const char* curr = p.path.data();
            const char* end = curr + p.path.size();

            while (curr < end) {
                const char* next_comma = curr;
                while (next_comma < end && *next_comma != ',') next_comma++;

                size_t vtx_id = 0;
                for (const char* ptr = curr; ptr < next_comma - 1; ++ptr) {
                    vtx_id = vtx_id * 10 + (*ptr - '0');
                }

                starts[out_idx] = current_pos + shift;
                orientedIds[out_idx] = vtx_id*((*(next_comma - 1) == '+') ? 1 : -1);

                current_pos += vtx_lengths[vtx_id];
                out_idx++;
                curr = next_comma + 1;
            }
        }

        // Step 4
        size_t next_vtx_id = vtx_lengths.size();
        for (size_t i = 0; i < data.idToSeq.size(); ++i) {
            if (counts[i] == 0) {
                size_t out_idx = offsets[i];
                starts[out_idx] = data.startingPos[i];
                orientedIds[out_idx] = next_vtx_id++;
            }
        }

        if (!data.startingPos.empty()) {
            starts.push_back(data.startingPos.back());
        }
    }
};

struct GraphB {
    std::vector<std::vector<long>> occurrences;
    std::vector<size_t> lengths, labelStartingPos;
    std::string labels = "";

    void init(const std::string& gfa, FastaData& data) {
        std::ifstream file(gfa);
        if (!file.is_open()) return;

        std::string line;
        std::vector<long> starts, ids;
        lengths.push_back(0), labelStartingPos.push_back(0);
        long labelPos = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == 'S') {
                size_t t1 = line.find('\t');
                size_t t2 = line.find('\t', t1 + 1);
                size_t t3 = line.find('\t', t2 + 1);

                size_t seq_end = (t3 == std::string::npos) ? line.length() : t3;
                std::string seq = line.substr(t2 + 1, seq_end - (t2 + 1));

                labelStartingPos.push_back(labelPos);
                lengths.push_back(seq.length());
                labels.append(seq);
                labelPos += seq.length();

            } else if (line[0] == 'P') {
                size_t t1 = line.find('\t');
                size_t t2 = line.find('\t', t1 + 1);
                size_t t3 = line.find('\t', t2 + 1);

                if (t2 == std::string::npos) continue; 
                if (t3 == std::string::npos) t3 = line.length();

                std::string genome = line.substr(t1 + 1, t2 - t1 - 1);

                auto it = data.seqToId.find(genome);
                if (it == data.seqToId.end()) continue;
                long shift = data.startingPos[it->second];

                const char* path_start = line.data() + t2 + 1;
                const char* path_end = line.data() + t3;
                const char* current = path_start;
                long current_pos = 0;

                while (current < path_end) {
                    const char* next_comma = current;
                    while (next_comma < path_end && *next_comma != ',') {
                        next_comma++;
                    }

                    if (next_comma - current < 2) break; // Zabezpieczenie przed pustym vtx

                    long vtx_id = 0;

                    for (const char* p = current; p < next_comma - 1; ++p) {
                        vtx_id = vtx_id * 10 + (*p - '0');
                    }

                    int orient = (*(next_comma - 1) == '+') ? 1 : -1;

                    starts.push_back((current_pos + shift) * orient);
                    ids.push_back(vtx_id);
                    current_pos += lengths[vtx_id];

                    current = next_comma + 1;
                }
            }
        }

        occurrences.assign(lengths.size(), std::vector<long>());
        for (size_t i = 0; i < ids.size(); i++) {
            occurrences[ids[i]].push_back(starts[i]);
        }
    }
};

int main(int argc, char* argv[]) {
        FastaData data;
        GraphA graphA;
        GraphB graphB;
        data.init(argv[3]), graphB.init(argv[2], data), graphA.init(argv[1], data);
        std::ofstream out(argv[4]);
        std::cout << data.startingPos.size() << " " << data.idToSeq.size() << "\n";
        std::cout << graphB.occurrences.size() << " " << graphB.lengths.size() << "\n";
        std::set<int> breakpoints;
        std::map<std::pair<size_t, int>, int> match, newVtxId;
        std::vector<std::pair<size_t, int>> matchingVector;
        std::vector<int> g_ids;
        size_t vtxId = 1;

        for (int i = 1; i < graphB.occurrences.size(); i++) {
                if (graphB.occurrences[i].size()<=1) continue;
                breakpoints.clear();
                breakpoints.insert(0), breakpoints.insert(graphB.lengths[i]);
                g_ids.clear();

                for (const auto& start : graphB.occurrences[i]) {
                        auto g_it = std::lower_bound(data.startingPos.begin(), data.startingPos.end(), std::abs(start));
                        if (*g_it != std::abs(start)) --g_it;
                        int g_id = std::distance(data.startingPos.begin(), g_it);
                        g_ids.push_back(g_id);
                        auto it = std::lower_bound(graphA.starts.begin(), graphA.starts.end(), std::abs(start));
                        size_t start_idx = std::distance(graphA.starts.begin(), it);

                        it = std::lower_bound(graphA.starts.begin(), graphA.starts.end(), std::abs(start)+graphB.lengths[i]);
                        size_t end_idx = std::distance(graphA.starts.begin(), it);

                        for (int k = start_idx; k < end_idx; k++) {
                                int strand = (start>0) ? 1 : -1;
                                if (strand==1) breakpoints.insert(graphA.starts[k]-start);
                                else {breakpoints.insert(graphB.lengths[i]-start-graphA.starts[k]);}
                        }

                }

                for (auto it=breakpoints.begin(); std::next(it)!=breakpoints.end(); ++it) {
                        auto currentBp = *it, nextBp = *std::next(it);
                        int len = nextBp-currentBp;
                        match.clear(), matchingVector.clear(), newVtxId.clear();
                        std::string subLabel = graphB.labels.substr(graphB.labelStartingPos[i]+currentBp, len);

                        for (const auto& start : graphB.occurrences[i]) {
                                auto shift = (start>0) ? currentBp : graphB.lengths[i]-nextBp;
                                auto itt = std::lower_bound(graphA.starts.begin(), graphA.starts.end(), std::abs(start)+shift);
                                if (*itt != std::abs(start)+shift) --itt;
                                size_t start_idx = std::distance(graphA.starts.begin(), itt);

                                if (graphA.orientedIds[start_idx]>0) {
                                        match[{graphA.orientedIds[start_idx], std::abs(start)+shift-graphA.starts[start_idx]}]++;
                                        matchingVector.push_back({graphA.orientedIds[start_idx], std::abs(start)+shift-graphA.starts[start_idx]});
                                } else {
                                        match[{(-1)*graphA.orientedIds[start_idx], graphA.starts[start_idx+1] - std::abs(start)-shift-len}]++;
                                        matchingVector.push_back({(-1)*graphA.orientedIds[start_idx], graphA.starts[start_idx+1] - std::abs(start)-shift-len});
                                }
                        }
                        for (const auto& [key, val]: match) {
                                if (val<=1) continue;
                                newVtxId[key] = vtxId;
                                out << "S\t" << vtxId << "\t" << subLabel << "\n";
                                vtxId++;
                        }

                        for (int j = 0; j < matchingVector.size(); j++) {
                                if (match[matchingVector[j]]<=1) continue;
                                char strand = (graphB.occurrences[i][j]>0) ? '+' : '-';
                                size_t from = std::abs(graphB.occurrences[i][j])-data.startingPos[g_ids[j]]+currentBp, to = from+len;
                                out << "P\t" << data.idToSeq[g_ids[j]] << ":" << from << "-" << to << "\t" << newVtxId[matchingVector[j]] << strand << "\n"; 
                        }
                }
        }
        out.close();
        return 0;
}

