#ifndef _DSTREAM_DSTREAM_H_
#define _DSTREAM_DSTREAM_H_

#include <stdint.h>
#include <unordered_map>
#include <vector>

#include "CharacteristicVector.h"
#include "Cluster.h"
#include "Key.h"

typedef std::unordered_map<Key, CharacteristicVector * > Gridlist;
typedef std::unordered_map<unsigned int, Cluster *> Clusters;

void AdjustClustering(Gridlist & grid_list, Clusters & clusters, uint64_t time_now, float d_m, float d_l);
void CalculateDensityParams(float & d_m, float & d_l);
void CallClusteringOnGrid(Key grid, Gridlist & grid_list, Clusters & clusters);
void CalculateXYCoords(double dx, double dy, double & x, double & y);
void CalculateXYDistances(double x, double y, double & dx, double & dy);
bool CheckForDenseNeighbors(Key & key, Gridlist & grid_list);
bool CheckIfConnected(Clusters & clusters, Gridlist & grid_list, uint32_t label);
void CheckNeighborsConnectivity(Key grid, Gridlist & grid_list, Clusters & clusters);
std::string CreateClustersJsonString(Clusters & clusters, Gridlist & grid_list);
std::string CreateCoordsJsonString(float lat, float lng);
std::string CreatePathJsonString(std::vector<Key> & path, uint32_t label);
std::string CreateGridsJsonString(Gridlist & grid_list);
float DensityThresholdFunction(uint64_t time_updated, uint64_t time_now);
void Deserialize(uint8_t * buffer, uint32_t buffer_size, uint64_t & time_now, Gridlist & grid_list);
float EstimatedDensitiesSum(uint64_t time_now);
void GenerateRandomPair(float & x, float & y, int counter, std::vector<Key> & points);
uint32_t GetAmountOfGridsInClusters(Clusters & clusters);
void GetNeighbors(float x, float y, Key neighbors[4]);
bool IsSporadic(float density, uint64_t time_updated, uint64_t time_now);
uint32_t LabelHash(Key & key);
void MergeChangesToBuffer(uint8_t * buffer, uint32_t buffer_size, Gridlist & grid_list);
void MergeClusters(Cluster & main_cluster, Cluster & secondary_cluster, Gridlist & grid_list, Clusters & clusters);
void PrintClusters(Clusters & clusters);
void PrintGridsByClusters(Clusters & clusters, Gridlist & grid_list);
void PrintPath(std::vector<Key> & path);
void PrintTable(Gridlist & grid_list, uint64_t time_now);
void ReassembleClusters(Gridlist & grid_list, Clusters & clusters);
void RelabelNeighbors(Key key, Gridlist & grid_list, uint32_t old_label, uint32_t new_label);
void RemoveEmptyClusters(Clusters & clusters);
void RemoveSporadic(Gridlist & grid_list, uint64_t time_now);
bool RenameCluster(Clusters & clusters, Gridlist & grid_list, uint32_t label);
uint8_t * SerializeClusters(Clusters & clusters, Gridlist & grid_list, uint32_t & buffer_size);
uint8_t * SerializeGrids(Gridlist & grid_list, uint32_t & buffer_size);
void SplitCluster(Clusters & clusters, Gridlist & grid_list, std::vector<int> checklist, uint32_t label, int count);
void UpdateDensities(Gridlist & grid_list, uint64_t time_now, float d_m, float d_l);

#endif // !_DSTREAM_DSTREAM_H_