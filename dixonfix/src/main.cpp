/*
* dixonfix - command line tool for correcting fat-water swaps in Dixon MRI.
*
* Copyright 2016 Ben Glocker <b.glocker@imperial.ac.uk>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
* Please cite:
*  B. Glocker, E. Konukoglu, I. Lavdas, J.E. Iglesias, E.O. Aboagye, A.G. Rockall, D. Rueckert
*  Correction of Fat-Water Swaps in Dixon MRI
*  International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 2016
*
*/

#include "RandomField.h"
#include "Optimization.h"
#include "miaImage.h"
#include "miaImageProcessing.h"
#include "itkio.h"
#include "itkBilateralImageFilter.h"

#include <vector>
#include <iostream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace mrfopt;
using namespace mia;

class DixonFixEnergyFunction : public EnergyFunction<double, int>
{
public:
  DixonFixEnergyFunction(const Image &likelihoodMap0, const Image &likelihoodMap1, const std::vector<int> &nodeImageIndices, double lambda)
    : m_likelihood_0(likelihoodMap0)
    , m_likelihood_1(likelihoodMap1)
    , m_node_img_indices(nodeImageIndices)
    , m_lambda(lambda)
  {}

  double unary_potential(int nodeIndex, int label) const override
  {
    int img_idx = m_node_img_indices[nodeIndex];
    return (label == 0) ? m_likelihood_0.data()[img_idx] : m_likelihood_1.data()[img_idx];
  }

  double clique_potential(const Clique& clique, const std::vector<int>& labels) const override
  {
    auto clique_size = static_cast<int>(clique.size());
    if (clique_size == 2)
    {
      return pairwise_potential(clique.nodes[0], clique.nodes[1], labels[0], labels[1]);
    }
    else
    {
      return 0;
    }
  }

  double pairwise_potential(int nodeA, int nodeB, int labelA, int labelB) const
  {
    return (labelA == labelB) ? 0.0 : m_lambda;
  }

private:
  Image m_likelihood_0;
  Image m_likelihood_1;
  std::vector<int> m_node_img_indices;
  double m_lambda;
};

void likelihood_maps(const Image &f, const Image &w, const Image &fp, const Image &wp, const Image &fu, const Image &wu, const Image &labelmap, Image &map0, Image &map1)
{
  for (size_t i = 0; i < f.size(); i++)
  {
    if (labelmap.data()[i] > 0)
    {
      auto fu_i = std::max(fu.data()[i], 1.0f);
      auto wu_i = std::max(wu.data()[i], 1.0f);

      auto sqDiffF0 = abs(f.data()[i] - fp.data()[i]) / fu_i;
      auto sqDiffW0 = abs(w.data()[i] - wp.data()[i]) / wu_i;
      map0.data()[i] = sqDiffW0 + sqDiffF0;

      auto sqDiffF1 = abs(f.data()[i] - wp.data()[i]) / wu_i;
      auto sqDiffW1 = abs(w.data()[i] - fp.data()[i]) / fu_i;
      map1.data()[i] = sqDiffW1 + sqDiffF1;

      //map0.data()[i] = sqDiffW0 / (sqDiffW0 + sqDiffW1) + sqDiffF0 / (sqDiffF0 + sqDiffF1);
      //map1.data()[i] = sqDiffW1 / (sqDiffW0 + sqDiffW1) + sqDiffF1 / (sqDiffF0 + sqDiffF1);
    }
  }
}

void bilateral_filter(const Image &input, Image &output, double domainSigma, double rangeSigma)
{
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::BilateralImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  ImageType::Pointer itkInputImage = itkio::miaToItk(input);
  double domainSigmas[Dimension];
  domainSigmas[0] = domainSigma;
  domainSigmas[1] = domainSigma;
  domainSigmas[2] = domainSigma;
  filter->SetInput(itkInputImage);
  filter->SetDomainSigma(domainSigmas);
  filter->SetRangeSigma(rangeSigma);
  filter->Update();

  ImageType::Pointer itkOutputImage = filter->GetOutput();

  ImageType::RegionType region = itkOutputImage->GetLargestPossibleRegion();
  ImageType::IndexType startIndex;
  startIndex.Fill(0);
  region.SetIndex(startIndex);

  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  ConstIteratorType outputImageIt(itkOutputImage, region);
  size_t index = 0;
  while (!outputImageIt.IsAtEnd())
  {
    output.data()[index] = outputImageIt.Get();
    ++outputImageIt;
    index++;
  }
}

Graph masked_grid_graph_first_order(const Image &mask, std::vector<int> &nodeImageIndices)
{
  int nodesX = mask.sizeX();
  int nodesY = mask.sizeY();
  int nodesZ = mask.sizeZ();

  int num_nodes = 0;
  Image idx = mask.clone();
  fill(idx, -1);
  for (size_t i = 0; i < mask.size(); i++)
  {
    if (mask.data()[i] > 0)
    {
      idx.data()[i] = num_nodes;
      nodeImageIndices.push_back(i);
      num_nodes++;
    }
  }

  Graph graph;
  graph.num_nodes(num_nodes);

  for (int z = 0; z < nodesZ; z++)
  {
    for (int y = 0; y < nodesY; y++)
    {
      for (int x = 0; x < nodesX; x++)
      {
        int index = x + y * nodesX + z * nodesX * nodesY;
        if (mask.data()[index] > 0)
        {
          if (x < nodesX - 1)
          {
            if (mask.data()[index + 1] > 0)
            {
              graph.add_clique(Clique(idx.data()[index], idx.data()[index + 1]));
            }
          }
          if (y < nodesY - 1)
          {
            if (mask.data()[index + nodesX] > 0)
            {
              graph.add_clique(Clique(idx.data()[index], idx.data()[index + nodesX]));
            }
          }
          if (z < nodesZ - 1)
          {
            if (mask.data()[index + nodesX * nodesY] > 0)
            {
              graph.add_clique(Clique(idx.data()[index], idx.data()[index + nodesX* nodesY]));
            }
          }
        }
      }
    }
  }

  return graph;
}

void swap_labeling(const Image &f, const Image &w, const Image &fp, const Image &wp, const Image &fu, const Image &wu, Image &labelmap, double lambda, std::string folderOut, bool writeTemp)
{
  std::cout << "+ computing likelihood maps...";
  auto start = std::chrono::high_resolution_clock::now();

  Image map0 = f.clone();
  Image map1 = f.clone();
  zeros(map0);
  zeros(map1);
  likelihood_maps(f, w, fp, wp, fu, wu, labelmap, map0, map1);

  if (writeTemp)
  {
    map0.dataType(mia::FLOAT);
    map1.dataType(mia::FLOAT);

    std::stringstream filename_likelihood_0;
    filename_likelihood_0 << folderOut << "/likelihood_0.nii.gz";
    itkio::save(map0, filename_likelihood_0.str());
    std::stringstream filename_likelihood_1;
    filename_likelihood_1 << folderOut << "/likelihood_1.nii.gz";
    itkio::save(map1, filename_likelihood_1.str());
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "done. (" << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() << " ms)" << std::endl;

  std::cout << "+ running max-flow..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  std::vector<int> nodeImageIndices;
  Graph graph = masked_grid_graph_first_order(labelmap, nodeImageIndices);
  DixonFixEnergyFunction function(map0, map1, nodeImageIndices, lambda);

  std::vector<int> labeling(graph.num_nodes());
  std::fill(labeling.begin(), labeling.end(), 0);

  std::vector<int> proposal(graph.num_nodes());
  std::fill(proposal.begin(), proposal.end(), 1);

  auto energy_initial = compute_energy(graph, function, labeling);
  std::cout << " - initial energy:\t" << energy_initial << std::endl;

  int numImprovements = 0;
  Reduction reductionMode = HOCR;
  fusion_move(graph, function, labeling, proposal, numImprovements, reductionMode);

  auto energy_final = compute_energy(graph, function, labeling);
  std::cout << " - final energy:\t" << energy_final << std::endl;

  // generate swap labeling map
  int node_idx = 0;
  for (int i = 0; i < labelmap.size(); i++)
  {
    if (labelmap.data()[i] > 0)
    {
      labelmap.data()[i] = (labeling[node_idx] == 1) ? 1.0f : 0.0f;
      node_idx++;
    }
  }

  stop = std::chrono::high_resolution_clock::now();
  std::cout << " - finished. (" << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() << " ms)" << std::endl;
}

int main(int argc, char* argv[])
{
  std::string filename_F;
  std::string filename_W;
  std::string filename_F_pred;
  std::string filename_W_pred;
  std::string filename_F_std;
  std::string filename_W_std;
  std::string folder_output;
  double spacing;
  double lambda;
  //double domain_sigma;
  //double range_sigma;

  bool write_temp;

  try
  {
    // Declare the supported options.
    po::options_description options("options");
    options.add_options()
    ("help,h", "produce help message")
    ("temp,t", po::bool_switch(&write_temp)->default_value(false), "enable output of temporary files")
    ("fat,f", po::value<std::string>(&filename_F), "filename of fat image")
    ("water,w", po::value<std::string>(&filename_W), "filename of water image")
    ("fat_prediction,p", po::value<std::string>(&filename_F_pred), "filename of predicted fat image")
    ("water_prediction,q", po::value<std::string>(&filename_W_pred), "filename of predicted water image")
    ("fat_uncertainty,u", po::value<std::string>(&filename_F_std), "filename of predicted fat image uncertainty map")
    ("water_uncertainty,v", po::value<std::string>(&filename_W_std), "filename of predicted water image uncertainty map")
    ("output,o", po::value<std::string>(&folder_output), "foldername for output")
    ("spacing,s", po::value<double>(&spacing), "voxel spacing in mm used for computations")
    ("lambda,l", po::value<double>(&lambda), "lambda for spatial regularization")
    //("domain_sigma,d", po::value<double>(&domain_sigma), "domain sigma for bilateral filter")
    //("range_sigma,r", po::value<double>(&range_sigma), "range sigma for bilateral filter")
    ;

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help") || vm.size() == 0)
    {
      std::cout << options << std::endl;
      return 0;
    }
  }
  catch (std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return 1;
  }

  std::cout << std::endl;
  std::cout << "            DIXON-FIX" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  std::cout << "+ loading and resampling images...";
  auto start = std::chrono::high_resolution_clock::now();


  Image F = itkio::load(filename_F);
  Image W = itkio::load(filename_W);
  Image FP = itkio::load(filename_F_pred);
  Image WP = itkio::load(filename_W_pred);
  Image FU = itkio::load(filename_F_std);
  Image WU = itkio::load(filename_W_std);

  // resample all into native resolution of F
  Image W_res = F.clone();
  resample(W, W_res, mia::LINEAR);

  Image FP_res = F.clone();
  Image WP_res = F.clone();
  resample(FP, FP_res, mia::LINEAR);
  resample(WP, WP_res, mia::LINEAR);

  Image FU_res = F.clone();
  Image WU_res = F.clone();
  resample(FU, FU_res, mia::LINEAR);
  resample(WU, WU_res, mia::LINEAR);

  // low resolution images for first stage
  Image F_low_resampled = resample(F, Eigen::Vector3d(spacing, spacing, spacing), mia::LINEAR);
  Image F_low_res = subimage(F_low_resampled, 1, 1, 1, F_low_resampled.sizeX() - 2, F_low_resampled.sizeY() - 2, F_low_resampled.sizeZ() - 2).clone();
  Image W_low_res = F_low_res.clone();
  resample(W, W_low_res, mia::LINEAR);

  Image FP_low_res = F_low_res.clone();
  Image WP_low_res = W_low_res.clone();
  resample(FP, FP_low_res, mia::LINEAR);
  resample(WP, WP_low_res, mia::LINEAR);

  Image FU_low_res = F_low_res.clone();
  Image WU_low_res = W_low_res.clone();
  resample(FU, FU_low_res, mia::LINEAR);
  resample(WU, WU_low_res, mia::LINEAR);

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "done. (" << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() << " ms)" << std::endl;


  //std::cout << "+ smoothing images and predictions...";
// auto start = std::chrono::high_resolution_clock::now();

// bilateral_filter(f, f_filtered, domainSigma, rangeSigma);
// bilateral_filter(w, w_filtered, domainSigma, rangeSigma);
// bilateral_filter(fp, fp_filtered, domainSigma, rangeSigma);
// bilateral_filter(wp, wp_filtered, domainSigma, rangeSigma);
  //gauss(F_low_res, F_low_res, domainSigma, domainSigma, domainSigma);
  //gauss(W_low_res, W_low_res, domainSigma, domainSigma, domainSigma);
  //gauss(FP_res, FP_res, domainSigma, domainSigma, domainSigma);
  //gauss(WP_res, WP_res, domainSigma, domainSigma, domainSigma);

  //if (writeTemp)
  //{
  //  std::stringstream filename_F_bilateral;
  //  filename_F_bilateral << folderOut << "/F_bilateral.nii.gz";
  //  itkio::save(f_filtered, filename_F_bilateral.str());

  //  std::stringstream filename_W_bilateral;
  //  filename_W_bilateral << folderOut << "/W_bilateral.nii.gz";
  //  itkio::save(w_filtered, filename_W_bilateral.str());

  //  std::stringstream filename_FP_bilateral;
  //  filename_FP_bilateral << folderOut << "/FP_bilateral.nii.gz";
  //  itkio::save(fp_filtered, filename_FP_bilateral.str());

  //  std::stringstream filename_WP_bilateral;
  //  filename_WP_bilateral << folderOut << "/WP_bilateral.nii.gz";
  //  itkio::save(wp_filtered, filename_WP_bilateral.str());
  //}

  //auto stop = std::chrono::high_resolution_clock::now();
  //std::cout << "done. (" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms)" << std::endl;


  std::cout << std::endl;
  std::cout << "+ FIRST STAGE - LOW RESOLUTION" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  // compute swap labeling for low resolution images
  Image labelmap_low_res = F_low_res.clone();
  mia::ones(labelmap_low_res);
  swap_labeling(F_low_res, W_low_res, FP_low_res, WP_low_res, FU_low_res, WU_low_res, labelmap_low_res, lambda, folder_output, false);

  //if (write_temp)
  //{
  //  for (int i = 0; i < labelmap_low_res.size(); i++)
  //  {
  //    if (labelmap_low_res.data()[i] > 0.0f)
  //    {
  //      float temp = F_low_res.data()[i];
  //      F_low_res.data()[i] = W_low_res.data()[i];
  //      W_low_res.data()[i] = temp;
  //    }
  //  }

  //  std::stringstream filename_F_low_res_fixed;
  //  filename_F_low_res_fixed << folder_output << "/F_low_res_fixed.nii.gz";
  //  itkio::save(F_low_res, filename_F_low_res_fixed.str());

  //  std::stringstream filename_W_low_res_fixed;
  //  filename_W_low_res_fixed << folder_output << "/W_low_res_fixed.nii.gz";
  //  itkio::save(W_low_res, filename_W_low_res_fixed.str());
  //}

  std::cout << std::endl;
  std::cout << "+ SECOND STAGE - FULL RESOLUTION" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  // one voxel dilation on swap labeling
  dilate_binary(labelmap_low_res.clone(), labelmap_low_res);

  // compute constrainted swap labeling on native resolution of F
  Image labelmap = F.clone();
  resample(labelmap_low_res, labelmap, mia::NEAREST);
  swap_labeling(F, W_res, FP_res, WP_res, FU_res, WU_res, labelmap, lambda, folder_output, write_temp);


  std::cout << "+ applying swap labeling...";
  start = std::chrono::high_resolution_clock::now();

  // apply swap labeling
  for (int i = 0; i < labelmap.size(); i++)
  {
    if (labelmap.data()[i] > 0.0f)
    {
      float temp = F.data()[i];
      F.data()[i] = W_res.data()[i];
      W_res.data()[i] = temp;
    }
  }

  stop = std::chrono::high_resolution_clock::now();
  std::cout << "done. (" << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() << " ms)" << std::endl;


  std::cout << "+ saving results...";
  start = std::chrono::high_resolution_clock::now();

  if (!fs::exists(folder_output)) fs::create_directories(folder_output);

  std::stringstream filename_f_fixed;
  filename_f_fixed << folder_output << "/f_fixed.nii.gz";
  itkio::save(F, filename_f_fixed.str());

  std::stringstream filename_w_fixed;
  filename_w_fixed << folder_output << "/w_fixed.nii.gz";
  itkio::save(W_res, filename_w_fixed.str());

  std::stringstream filename_labelmap;
  filename_labelmap << folder_output << "/labelmap.nii.gz";
  itkio::save(labelmap, filename_labelmap.str());

  stop = std::chrono::high_resolution_clock::now();
  std::cout << "done. (" << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() << " ms)" << std::endl;

  std::cout << "--------------------" << std::endl;
}
