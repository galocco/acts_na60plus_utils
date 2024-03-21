import ROOT

def plot_efficiency_vs_method(directory, root_file_path1, root_file_path2, title1, title2, tree_name):
    # Open the ROOT files
    file1 = ROOT.TFile.Open(root_file_path1)
    file2 = ROOT.TFile.Open(root_file_path2)
    # Get the trees from the files
    efficiency1 = file1.Get(tree_name)
    efficiency2 = file2.Get(tree_name)

    # Create a TCanvas and draw the TEfficiency histograms
    canvas = ROOT.TCanvas("canvas", "Efficiency Comparison", 1600, 1200)
    efficiency1.SetLineColor(ROOT.kRed)
    efficiency1.SetMarkerColor(ROOT.kRed)
    efficiency2.SetLineColor(ROOT.kBlue)
    efficiency2.SetMarkerColor(ROOT.kBlue)
    efficiency1.SetMarkerStyle(20)
    efficiency2.SetMarkerStyle(21)
    efficiency2.Draw()
    ROOT.gPad.Update()
    graph = efficiency2.GetPaintedGraph()
    if "eff" in tree_name:
        graph.SetMinimum(0.5)
        graph.SetMaximum(1.02)
    else:
        graph.SetMinimum(-0.02)
        graph.SetMaximum(1.02)
    ROOT.gPad.Update()
    efficiency1.Draw("same")

    # Add a legend
    legend = ROOT.TLegend(0.8, 0.8, 1, 1)
    legend.AddEntry(efficiency1, title1, "lep")
    legend.AddEntry(efficiency2, title2, "lep")
    legend.Draw()
    # Show the canvas
    canvas.Draw()
    canvas.Print(directory+"/"+tree_name+"_comparison.png")

    # Close the ROOT files
    file1.Close()
    file2.Close()

def plot_efficiency_vs_step(directory, root_file_path1, root_file_path2, root_file_path3, title1, title2, title3, tree_name):
    # Open the ROOT files
    file1 = ROOT.TFile.Open(directory+root_file_path1)
    file2 = ROOT.TFile.Open(directory+root_file_path2)
    file3 = ROOT.TFile.Open(directory+root_file_path3)
    
    # Get the trees from the files
    efficiency1 = file1.Get(tree_name)
    efficiency2 = file2.Get(tree_name)
    efficiency3 = file3.Get(tree_name)

    # Create a TCanvas and draw the TEfficiency histograms
    canvas = ROOT.TCanvas("canvas", "Efficiency Comparison", 1600, 1200)
    efficiency1.SetLineColor(ROOT.kRed)
    efficiency1.SetMarkerColor(ROOT.kRed)
    efficiency2.SetLineColor(ROOT.kBlue)
    efficiency2.SetMarkerColor(ROOT.kBlue)
    efficiency1.SetMarkerStyle(20)
    efficiency2.SetMarkerStyle(21)
    efficiency3.SetMarkerStyle(22)
    efficiency3.Draw()
    ROOT.gPad.Update()
    graph = efficiency3.GetPaintedGraph()
    if "eff" in tree_name:
        graph.SetMinimum(0.5)
        graph.SetMaximum(1.02)
    else:
        graph.SetMinimum(-0.02)
        graph.SetMaximum(1.02)
    ROOT.gPad.Update()
    efficiency2.Draw("same")
    efficiency1.Draw("same")

    # Add a legend
    legend = ROOT.TLegend(0.8, 0.8, 1, 1)
    legend.AddEntry(efficiency1, title1, "lep")
    legend.AddEntry(efficiency2, title2, "lep")
    legend.AddEntry(efficiency3, title3, "lep")
    legend.Draw()
    # Show the canvas
    canvas.Draw()
    canvas.Print(directory+tree_name+"_comparison.png")

    # Close the ROOT files
    file1.Close()
    file2.Close()
    file3.Close()





var_list = ["pT","phi","eta"]
directory_list = ["/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun1/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun50/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun300/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun800/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_deadZones/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing_deadZones/",
                  "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing/"]
for var in var_list:
    for directory in directory_list:
        # Example usage
        plot_efficiency_vs_step(directory,
                        "performance_seeding.root",
                        "performance_ckf.root",
                        "performance_ambi.root",
                        "Seeding",
                        "CKF",
                        "Ambi",
                        "trackeff_vs_"+var)

step_list = ["seeding","ckf","ambi"]
rate_list = ["trackeff","fakerate","duplicationRate"]
for step in step_list:
    for var in var_list:
        for rate in rate_list:
            if step == "seeding" and rate == "fakerate":
                continue

            plot_efficiency_vs_method("figures",
                            "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_deadZones/performance_"+step+".root",
                            "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing_deadZones/performance_"+step+".root",
                            "Truth",
                            "Standard",
                            rate+"_vs_"+var)


def printResults(mode, file):
    import uproot
    rootFile = uproot.open(file)
    print("Standard")
    print("eff_track:",rootFile["eff_tracks"].member("fElements")[0])
    print("fake_track:",rootFile["fakerate_tracks"].member("fElements")[0])
    print("dupl_track:",rootFile["duplicaterate_tracks"].member("fElements")[0])

    print("eff_particle:",rootFile["eff_particles"].member("fElements")[0])
    print("fake_particle:",rootFile["fakerate_particles"].member("fElements")[0])
    print("dupl_particle:",rootFile["duplicaterate_particles"].member("fElements")[0])

    print("tot:",rootFile["eff_tracks"].member("fElements")[0]+rootFile["fakerate_tracks"].member("fElements")[0])
    print("tot:",rootFile["eff_particles"].member("fElements")[0]+rootFile["fakerate_particles"].member("fElements")[0])
#printResults("standard","/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing/performance_ckf.root")
#printResults("truth","/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing/performance_ckf.root")