import uproot
import ROOT
import math
import os

def getResolutions(file_path = "output_standardSeeding_iterativeVertexing/tracksummary_ambi.root", output="histogram.root",directory="Resolutions"):
    tree_name = "tracksummary"

    # Open the ROOT file and get the TTree
    root_file = uproot.open(file_path)
    tree = root_file[tree_name]

    # Create a ROOT file to store the histograms
    output_file = ROOT.TFile(output, "RECREATE")

    # Get the branch
    px_tree = tree["t_px"]
    py_tree = tree["t_py"]
    pz_tree = tree["t_pz"]
    p_tree = tree["t_p"]
    d0_tree = tree["t_d0"]
    z0_tree = tree["t_z0"]
    qop_tree = tree["eQOP_fit"]
    phi_tree = tree["ePHI_fit"]
    theta_tree = tree["eTHETA_fit"]
    d0_fit_tree = tree["eLOC0_fit"]
    z0_fit_tree = tree["eLOC1_fit"]
    hit_tree = tree["nMajorityHits"]
    chi2_tree = tree["chi2Sum"]
    ndf_tree = tree["NDF"]

    # Get the data from the branch as a numpy array
    px_data = px_tree.array(library="np")
    py_data = py_tree.array(library="np")
    pz_data = pz_tree.array(library="np")
    p_data = p_tree.array(library="np")
    qop_data = qop_tree.array(library="np")
    theta_data = theta_tree.array(library="np")
    phi_data = phi_tree.array(library="np")
    hit_data = hit_tree.array(library="np")
    d0_data = d0_tree.array(library="np")
    z0_data = z0_tree.array(library="np")
    d0_fit_data = d0_fit_tree.array(library="np")
    z0_fit_data = z0_fit_tree.array(library="np")

    ndf_data = ndf_tree.array(library="np")
    chi2_data = chi2_tree.array(library="np")


    # Create a histogram for the branch
    res_px = ROOT.TH1F("res_px", "residuals p_{x};p_{x rec}-p_{x gen} (GeV/c);entries", 1000, -1,1)
    res_py = ROOT.TH1F("res_py", "residuals p_{y};p_{y rec}-p_{y gen} (GeV/c);entries", 1000, -1,1)
    res_pz = ROOT.TH1F("res_pz", "residuals p_{z};p_{z rec}-p_{z gen} (GeV/c);entries", 1000, -1,1)
    res_pt = ROOT.TH1F("res_pt", "residuals p_{T};p_{T rec}-p_{T gen} (GeV/c);entries", 1000, -1,1)
    res_p = ROOT.TH1F("res_p", "residuals p;p_{rec}-p_{gen} (GeV/c);entries", 1000, -1,1)
    res_d0 = ROOT.TH1F("res_d0", "residuals d_{0};d_{0 rec}-d_{0 gen} (mm);entries", 1000, -1,1)
    res_z0 = ROOT.TH1F("res_z0", "residuals z_{0};z_{0 rec}-z_{0 gen} (mm);entries", 1000, -5,5)
    res_px_vs_px = ROOT.TH2F("res_px_vs_px", "residuals p_{x};p_{x} (GeV/c);p_{x rec}-p_{x gen} (GeV/c);entries",1000,-3,3, 1000, -3,3)
    res_py_vs_py = ROOT.TH2F("res_py_vs_py", "residuals p_{y};p_{y} (GeV/c);p_{y rec}-p_{y gen} (GeV/c);entries",1000,-3,3, 1000, -3,3)
    res_dpx_vs_dpy = ROOT.TH2F("res_dpx_vs_dpy", "residuals p_{x} vs residuals p_{y};p_{y rec}-p_{y gen} (GeV/c);py_{rec}-py_{gen} (GeV/c);entries",200,-3,3, 200, -3,3)
    res_pz_vs_pz = ROOT.TH2F("res_pz_vs_pz", "residuals p_{z};pz (GeV/c);p_{z rec}-p_{zgen} (GeV/c);entries",200,0,20, 1000, -1,1)
    res_pt_vs_pt = ROOT.TH2F("res_pt_vs_pt", "residuals p_{t};pt (GeV/c);p_{t rec}-p_{tgen} (GeV/c);entries",30,0,3, 1000, -1,1)
    res_p_vs_p = ROOT.TH2F("res_p_vs_p", "residuals p;p (GeV/c);p_{rec}-p_{gen} (GeV/c);entries",200,0,20, 1000, -1,1)
    res_d0_vs_p = ROOT.TH2F("res_d0_vs_p", "residuals d_{0};p (GeV/c);d_{0 rec}-d_{0 gen} (mm);entries",100,0,100, 1000, -1,1)
    res_d0_vs_pt = ROOT.TH2F("res_d0_vs_pt", "residuals d_{0};p_{T} (GeV/c);d_{0 rec}-d_{0 gen} (mm);entries",30,0,3, 1000, -1,1)
    res_z0_vs_p = ROOT.TH2F("res_z0_vs_p", "residuals z_{0};p (GeV/c);z_{0 rec}-z_{0 gen} (mm);entries",100,0,100, 100, -5,5)
    res_z0_vs_pt = ROOT.TH2F("res_z0_vs_pt", "residuals z_{0};p_{T} (GeV/c);z_{0 rec}-z_{0 gen} (mm);entries",30,0,3, 100, -5,5)
    chi2ndof = ROOT.TH1F("chi2ndof", "reduced chi2;chi2/NDoF;entries", 1000, 0,100)
    chi2ndof = ROOT.TH1F("chi2ndof", "reduced chi2;chi2/NDoF;entries", 1000, 0,100)
    #px = r*math.sin(theta)*math.cos(phi)
    nev=0
    

    # Check if the directory already exists
    if not os.path.exists(directory):
        # Create the directory
        os.makedirs(directory)
        print(f"Directory '{directory}' created successfully.")
    else:
        print(f"Directory '{directory}' already exists.")

    cv = ROOT.TCanvas("cv","cv")
    for px_ev,py_ev,pz_ev,p_ev,qop_ev,theta_ev,phi_ev,hit_ev,d0_ev,z0_ev,d0_fit_ev,z0_fit_ev,chi2_ev,ndf_ev in zip(px_data,py_data,pz_data,p_data,qop_data,theta_data,phi_data,hit_data,d0_data,z0_data,d0_fit_data,z0_fit_data,chi2_data,ndf_data):
        nev += 1

        for px,py,pz,p,qop,theta,phi,hit,z0,d0,d0_fit,z0_fit,chi2,ndf in zip(px_ev,py_ev,pz_ev,p_ev,qop_ev,theta_ev,phi_ev,hit_ev,d0_ev,z0_ev,d0_fit_ev,z0_fit_ev,chi2_ev,ndf_ev):
            if hit >=0:
                r = abs(1./qop)
                px_fit = r*math.sin(theta)*math.cos(phi)
                py_fit = r*math.sin(theta)*math.sin(phi)
                pz_fit = r*math.cos(theta)
                pt_fit = math.sqrt(px_fit**2+py_fit**2)
                pt = math.sqrt(px**2+py**2)
                p_fit = math.sqrt(px_fit**2+py_fit**2+pz_fit**2)
                p = math.sqrt(px**2+py**2+pz**2)
                chi2ndof.Fill(chi2/ndf)
                res_px.Fill(px_fit-px)
                res_py.Fill(py_fit-py)
                res_pz.Fill(pz_fit-pz)
                res_pt.Fill(pt_fit-pt)
                res_d0.Fill(d0_fit-d0)
                res_z0.Fill(z0_fit-z0)
                res_p.Fill(p_fit-p)
                res_px_vs_px.Fill(px,px_fit-px)
                res_py_vs_py.Fill(py,py_fit-py)
                res_pz_vs_pz.Fill(pz,pz_fit-pz)
                res_pt_vs_pt.Fill(pt,pt_fit-pt)
                res_d0_vs_p.Fill(p,d0_fit-d0)
                res_z0_vs_p.Fill(p,z0_fit-z0)
                res_p_vs_p.Fill(p,p_fit-p)

                #res_py_vs_px.Fill(py_fit,px_fit)
                #res_px_vs_py.Fill(px_fit,py_fit)
                res_dpx_vs_dpy.Fill(px-px_fit,py-py_fit)

    chi2ndof.Write()
    res_px.Write()
    res_py.Write()
    res_pz.Write()
    res_pt.Write()
    res_d0.Write()
    res_z0.Write()
    res_p.Write()
    res_p_vs_p.Write()
    res_px_vs_px.Write()
    res_py_vs_py.Write()
    res_pz_vs_pz.Write()
    res_pt_vs_pt.Write()
    res_d0_vs_p.Write()
    res_z0_vs_p.Write()
    res_dpx_vs_dpy.Write()
    print(directory)
    res_px.Draw("colz")
    cv.SaveAs(directory+"/res_px.png")
    res_py.Draw("colz")
    cv.SaveAs(directory+"/res_py.png")
    res_pz.Draw("colz")
    cv.SaveAs(directory+"/res_pz.png")
    res_pt.Draw("colz")
    cv.SaveAs(directory+"/res_pt.png")
    res_d0.Draw("colz")
    cv.SaveAs(directory+"/res_d0.png")
    res_z0.Draw("colz")
    cv.SaveAs(directory+"/res_z0.png")
    res_p.Draw("colz")
    cv.SaveAs(directory+"/res_p.png")
    res_p_vs_p.Draw("colz")
    cv.SaveAs(directory+"/res_p_vs_p.png")
    res_px_vs_px.Draw("colz")
    cv.SaveAs(directory+"/res_px_vs_px.png")
    res_py_vs_py.Draw("colz")
    cv.SaveAs(directory+"/res_py_vs_py.png")
    res_pz_vs_pz.Draw("colz")
    cv.SaveAs(directory+"/res_pz_vs_pz.png")
    res_pt_vs_pt.Draw("colz")
    cv.SaveAs(directory+"/res_pt_vs_pt.png")
    res_d0_vs_p.Draw("colz")
    cv.SaveAs(directory+"/res_d0_vs_pt.png")
    res_z0_vs_p.Draw("colz")
    cv.SaveAs(directory+"/res_z0_vs_pt.png")
    res_dpx_vs_dpy.Draw("colz")
    cv.SaveAs(directory+"/res_dpx_vs_dpy.png")


    # Close the ROOT file and output file
    output_file.Close()
    root_file.close()

if False:
    directory_list = ["/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun1/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun50/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun300/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_gun800/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_truthEstimated_truthVertexing_deadZones/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing_deadZones/",
                    "/home/giacomo/ACTS_for_NA60+/acts_na60plus_utils/output_standardSeeding_iterativeVertexing/"]
    suffix_list = ["truth_gun1",
                    "truth_gun50",
                    "truth_gun300",
                    "truth_gun800",
                    "truth_deadZones",
                    "standard_deadZones",
                    "standard"]

    for directory, suffix in zip(directory_list, suffix_list):
        getResolutions(file_path = directory+"tracksummary_ambi.root", output="histogram"+suffix+".root",directory="Resolutions"+suffix)
    #getResolutions(file_path = "output_truthEstimated_truthVertexing/tracksummary_ambi.root", output="histogramTruth.root",directory="ResolutionsTruth")
    #getResolutions(file_path = "output_particleSmearing_truthVertexing/tracksummary_ambi.root", output="histogramSmeareds.root",directory="ResolutionsSmeareds")
    


def getD0Mass(file_path = "output_truthEstimated_truthVertexing/tracksummary_ckf.root",
              particle_path = "output_truthEstimated_truthVertexing/particles.root",
              output="histogramMass.root",
              directory="ResolutionsMass"):
    mD0 = 1.864840
    # Open the ROOT file and get the TTree
    particle_file = uproot.open(particle_path)
    particle_tree = particle_file["particles"]
    process_tree = particle_tree["process"]
    particles_id_tree = particle_tree["particle_id"]
    particles_type_tree = particle_tree["particle_type"]

    process_data = process_tree.array(library="pd")
    particles_id_data = particles_id_tree.array(library="pd")
    particles_type_data = particles_type_tree.array(library="pd")
    list_id = []
    for proc_ev, id_ev, type_ev in zip(process_data, particles_id_data, particles_type_data):
        pair = []
        for _ in range(0,10):
            pair.append([0,0])
        for proc, id, pdg in zip(proc_ev, id_ev, type_ev):
            #print(proc)
            if proc>0:
                if pdg == 211:
                    pair[proc][0] = id
                elif pdg == 321:
                    pair[proc][1] = id
        list_id.append(pair)
    # Open the ROOT file and get the TTree
    root_file = uproot.open(file_path)
    tree = root_file["tracksummary"]

    # Create a ROOT file to store the histograms
    output_file = ROOT.TFile(output, "RECREATE")

    # Get the branch
    qop_tree = tree["eQOP_fit"]
    phi_tree = tree["ePHI_fit"]
    theta_tree = tree["eTHETA_fit"]
    hit_tree = tree["nMajorityHits"]
    id_tree = tree["majorityParticleId"]

    # Get the data from the branch as a numpy array
    qop_data = qop_tree.array(library="np")
    theta_data = theta_tree.array(library="np")
    phi_data = phi_tree.array(library="np")
    hit_data = hit_tree.array(library="np")

    id_data = id_tree.array(library="np")


    # Create a histogram for the branch
    hist_mass = ROOT.TH1F("hist_mass", "D^{0} mass distribution;m_{#piK} (GeV/c);entries", 500, -1+mD0,1+mD0)
    hist_all_vs_pt = ROOT.TH2F("hist_all_vs_pt", "residuals mass;p_{T} GeV/#it{c}; m_{rec} (GeV/c);entries",15,0,3, 500, -1+mD0,1+mD0)
    hist_sigma_mass = ROOT.TH1F("hist_sigma_mass","D^{0} mass resolution;p_{T} GeV/#it{c};#sigma_{m} (MeV/c)",15,0,3)
    hist_std_mass = ROOT.TH1F("hist_std_mass","D^{0} mass resolution;p_{T} GeV/#it{c};RMS_{m} (MeV/c)",15,0,3)
    hist_all = ROOT.TH1F("hist_all", "residuals mass;m_{rec}-m_{pdg} (GeV/c);entries", 1000, -1+mD0,1+mD0)
    res_mass = ROOT.TH1F("res_mass", "residuals mass;m_{rec}-m_{pdg} (GeV/c);entries", 250, -0.5,0.5)
    res_mass_vs_pt = ROOT.TH2F("res_mass_vs_pt", "residuals mass;p_{T} GeV/#it{c}; m_{rec}-m_{pdg} (GeV/c);entries",15,0,3, 250, -0.5,0.5)
    res_mass_8 = ROOT.TH1F("res_mass_8", "residuals mass;m_{rec}-m_{pdg} (GeV/c);entries", 1000, -0.5,0.5)
    res_mass_9 = ROOT.TH1F("res_mass_9", "residuals mass;m_{rec}-m_{pdg} (GeV/c);entries", 1000, -0.5,0.5)
    res_mass_10 = ROOT.TH1F("res_mass_10", "residuals mass;m_{rec}-m_{pdg} (GeV/c);entries", 1000, -0.5,0.5)

    nev=0
    # Check if the directory already exists
    if not os.path.exists(directory):
        # Create the directory
        os.makedirs(directory)
        print(f"Directory '{directory}' created successfully.")
    else:
        print(f"Directory '{directory}' already exists.")

    cv = ROOT.TCanvas("cv","cv",1600,1200)
    v1 = ROOT.TLorentzVector(0,0,0,0)
    v2 = ROOT.TLorentzVector(0,0,0,0)
    for qop_ev,theta_ev,phi_ev,hit_ev,id_ev in zip(qop_data,theta_data,phi_data,hit_data,id_data):
        print("event: ",nev)
        if nev == -1:
            break
        nhit=0
        vmom_list = []
        cont_list = []
        for _ in range(0,10):
            vmom_list.append(ROOT.TLorentzVector(0,0,0,0))
            cont_list.append(0)

        for qop,theta,phi,hit,id in zip(qop_ev,theta_ev,phi_ev,hit_ev,id_ev):
            v = ROOT.TLorentzVector()
            index = 0
            ch1 = 1 if qop>0 else -1
 
            if False:
                r = abs(1./qop)
                px = r*math.sin(theta)*math.cos(phi)
                py = r*math.sin(theta)*math.sin(phi)
                pz = r*math.cos(theta)
                v1.SetXYZM(px,py,pz,0.139570)
                for qop2,theta2,phi2 in zip(qop_ev,theta_ev,phi_ev):
                    ch2 = 1 if qop2>0 else -1
                    if ch2*ch1 > 0:
                        continue

                    r2 = abs(1./qop2)
                    px2 = r2*math.sin(theta2)*math.cos(phi2)
                    py2 = r2*math.sin(theta2)*math.sin(phi2)
                    pz2 = r2*math.cos(theta2)
                    v2.SetXYZM(px2,py2,pz2,0.493677); 
                    v2 += v1
                    hist_all.Fill(v2.M())
                    hist_all_vs_pt.Fill(v2.Pt(),v2.M())
            for i in range(0,10):
                if id==list_id[nev][i][0]:
                    mass = 0.139570
                    index = i
                elif id == list_id[nev][i][1]:
                    mass = 0.493677
                    index = i
            if index != 0:
                nhit+=hit
                r = abs(1./qop)
                px = r*math.sin(theta)*math.cos(phi)
                py = r*math.sin(theta)*math.sin(phi)
                pz = r*math.cos(theta)
                v.SetXYZM(px,py,pz,mass); 
                vmom_list[index] += v
                cont_list[index] += 1
        for vmom,cont in zip(vmom_list,cont_list):
            if cont==2:
                res_mass.Fill(vmom.M()-mD0)
                hist_mass.Fill(vmom.M())
                
                res_mass_vs_pt.Fill(vmom.Pt(),vmom.M()-mD0)
                if nhit==8:
                    res_mass_8.Fill(vmom.M()-mD0)
                if nhit==9:
                    res_mass_9.Fill(vmom.M()-mD0)
                if nhit==10:
                    res_mass_10.Fill(vmom.M()-mD0)
        nev += 1

    res_mass.Draw("colz")
    cv.SaveAs("mass.png")
    res_mass.Write()
    hist_mass.Write()
    res_mass_vs_pt.Write()
    res_mass_8.Write()
    res_mass_9.Write()
    res_mass_10.Write()
    hist_all.Write()
    hist_all_vs_pt.Write()

    tf1 = ROOT.TF1("gaus","gausn(0)")
    for iBin in range(1,res_mass_vs_pt.GetNbinsX()+1):
        hist = res_mass_vs_pt.ProjectionY("hist",iBin,iBin)
        tf1.SetParameter(0,hist.GetEntries()*hist.GetBinWidth(2))
        tf1.SetParLimits(0,hist.GetEntries()*hist.GetBinWidth(2)/4,hist.GetEntries()*hist.GetBinWidth(2)*10)
        tf1.SetParameter(1,0)
        tf1.SetParameter(2,0.015)
        tf1.SetNpx(400)
        hist.Fit(tf1,"MQ+")
        hist.Draw()
        hist.SetTitle(str(round(res_mass_vs_pt.GetXaxis().GetBinLowEdge(iBin),1))+" #leq #it{p}_{T} < "+str(round(res_mass_vs_pt.GetXaxis().GetBinUpEdge(iBin),1))+" GeV/#it{c}")
        cv.SaveAs("test"+str(iBin)+".png")
        hist_std_mass.SetBinContent(iBin,hist.GetRMS()*1000)
        hist_sigma_mass.SetBinContent(iBin,tf1.GetParameter(2)*1000)
        hist_std_mass.SetBinError(iBin,hist.GetRMSError()*1000)
        hist_sigma_mass.SetBinError(iBin,tf1.GetParError(2)*1000)
    hist_std_mass.Write()
    hist_sigma_mass.Write()

    hist_std_mass.SetLineColor(ROOT.kRed)
    hist_std_mass.SetMarkerColor(ROOT.kRed)
    hist_std_mass.SetMarkerStyle(20)
    hist_sigma_mass.SetMarkerStyle(21)
    hist_std_mass.GetYaxis().SetTitle("Mass resolution (MeV/#it{c})")
    hist_std_mass.Draw()
    hist_sigma_mass.Draw("same")
    hist_std_mass.GetYaxis().SetRangeUser(0,80)
    legend = ROOT.TLegend(0.5,0.75,0.7,0.9)
    legend.AddEntry(hist_std_mass, "RMS", "lep")
    legend.AddEntry(hist_sigma_mass, "#sigma", "lep")
    legend.Draw()
    cv.SaveAs("D0resolution.png")
    

    # Close the ROOT file and output file
    output_file.Close()
    root_file.close()

getD0Mass()