import ROOT
import os
ROOT.gROOT.SetBatch(True)


def save_root_objects_as_png(directory, root_filename, output_dir):
    # Open the ROOT file
    root_file = ROOT.TFile.Open(directory+root_filename)
    if not root_file:
        print(f"Error: Failed to open ROOT file '{root_filename}'")
        return
    output_dir = directory+output_dir
    # Check if the directory already exists
    if not os.path.exists(output_dir):
        # Create the directory
        os.makedirs(output_dir)
        print(f"Directory '{output_dir}' created successfully.")
    else:
        print(f"Directory '{output_dir}' already exists.")

    # Loop over all keys in the file
    keys = [key.GetName() for key in root_file.GetListOfKeys()]
    for key_name in keys:
        obj = root_file.Get(key_name)
        if not obj:
            print(f"Error: Failed to retrieve object with key '{key_name}'")
            continue

        # Convert the object to an image and save it as PNG
        canvas = ROOT.TCanvas("canvas", "Canvas", 1600, 1200)
        obj.Draw()
        canvas.Update()

        image_name = f"{output_dir}/{key_name}.png"
        canvas.SaveAs(image_name)

        # Close the canvas to avoid memory leaks
        canvas.Close()

        print(f"Saved object '{key_name}' as {image_name}")

    # Close the ROOT file
    root_file.Close()

# Example usage:
root_filename = "performance_track_fitter_ckf.root"
output_directory = "figurePulls"


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
    save_root_objects_as_png(directory,root_filename, output_directory)
