import pandas as pd
import os
import sys
import math


class track_record:
    def __init__(self, start_frame, track):
        self.start_frame = start_frame
        self.track = track

    def add_cell(self, cell_id):
        self.track.append(cell_id)

    def get_last_cell(self):
        return self.track[len(self.track) - 1]
    
    def get_valid_track_len(self):
        return len([i for i in self.track if i > 0])

if __name__ == '__main__':
    base_dir = sys.argv[1]
    cell_sum_file = os.path.join(base_dir,'cell_summary.csv')
    df_cell_summaries = pd.read_csv(cell_sum_file, sep=',')
    frame_count = df_cell_summaries["z_id"].max()
    track_lists = []
    for i in range(1,frame_count+1):
        df_this_frame_cells = df_cell_summaries.loc[df_cell_summaries['z_id']==i].copy()
        pending_cells = df_this_frame_cells['cell_id'].tolist()
        min_assigned_cells = {}
        for j, track_rec in enumerate(track_lists):
            last_cell_id = track_rec.get_last_cell()
            if last_cell_id != 0:
                last_cell_cent_x = df_cell_summaries.loc[(df_cell_summaries['z_id']==i-1) & (df_cell_summaries['cell_id']==last_cell_id)]['centroid_x'].iloc[0]
                last_cell_cent_y = df_cell_summaries.loc[(df_cell_summaries['z_id']==i-1) & (df_cell_summaries['cell_id']==last_cell_id)]['centroid_y'].iloc[0]
                min_cell_id = 5000
                min_cell_dist = 5000
                for index, row in df_this_frame_cells.iterrows():
                    cell_id = int(row["cell_id"])
                    cent_x = row["centroid_x"]
                    cent_y = row["centroid_y"]
                    euclidean_distance = math.sqrt((last_cell_cent_x - cent_x) ** 2 + (last_cell_cent_y - cent_y) ** 2)
                    if euclidean_distance < min_cell_dist:
                        min_cell_dist=euclidean_distance
                        min_cell_id=cell_id

                if min_cell_dist > 100:
                    track_rec.add_cell(0)
                else:
                    min_assigned_cells[j] = [min_cell_id, min_cell_dist]
                    if min_cell_id in pending_cells:
                        pending_cells.remove(min_cell_id)
            else:
                track_rec.add_cell(0)

        # Assign cell to track with min dist
        for j, track_rec in enumerate(track_lists):
            if j in min_assigned_cells:
                next_cell = min_assigned_cells[j][0]
                next_cell_dist = min_assigned_cells[j][1]
                found_alt_min=False
                for key, value in min_assigned_cells.items():
                    if value[0]==next_cell and value[1]<next_cell_dist:
                        found_alt_min=True
                        break
                if not found_alt_min:
                    track_rec.add_cell(next_cell)
                else:
                    track_rec.add_cell(0)
        ##Start the new ones not assigned to any current tracks
        for cell_id in pending_cells:
            track_rec = track_record(i,[cell_id])
            track_lists.append(track_rec)

    min_tracks = 3
    output_file = os.path.join(base_dir,'cell_tracks.csv')
    with open(output_file, "w", newline="") as f:
        for i in range(1, frame_count + 1):
            f.write("Frame_" + str(i)+",")
        f.write("\n")
        for track_rec in track_lists:
            if track_rec.get_valid_track_len() > min_tracks:
                for i in range(1,track_rec.start_frame):
                    f.write("0,")
                track_str = ','.join([str(i) for i in track_rec.track])
                f.write(track_str)
                f.write("\n")