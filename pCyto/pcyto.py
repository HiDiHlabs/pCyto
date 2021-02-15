import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import matplotlib.patches as patches

giemsa_color_map = {"gpos100": (0., 0., 0.),
					"gpos": (0., 0., 0.),
					"gpos75": (130./255., 130./255., 130./255.),
					"gpos66": (160./255., 160./255., 160./255.),
					"gpos50": (200./255., 200./255., 200./255.),
					"gpos33": (210./255., 210./255., 210./255.),
					"gpos25": (200./255., 200./255., 200./255.),
					"gvar": (220./255., 220./255., 220./255.),
					"gneg": (255./255., 255./255., 255./255.),
					"acen": (217./255., 47./255., 39./255.),
					"stalk": (100./255., 127./255., 164./255.)}

def bezierMidPoint(x, x0, x2, t):
	return (((1-t)**2)*x0+(t**2)*x2-x)/(-2*(1-t)*t)

def plotFragments(fragment_list, clip_patch, color_map, alignment="vertical", ax = None):
	ax = ax if ax is not None else plt.gca

	for fragment in fragment_list:
		start = fragment[1]
		end = fragment[2]
		giemsa_stain = fragment[4]
		rect = None
		if(alignment == "vertical"):
			rect = Rectangle([0, -1*end], 1, end-start, facecolor=color_map[giemsa_stain], edgecolor="none", capstyle="butt")
		else:
			rect = Rectangle([start, 0], end-start, 1, facecolor=color_map[giemsa_stain], edgecolor="none", capstyle="butt")
		rect.set_clip_path(clip_patch)

		ax.add_patch(rect)


def plotChromosomeIdeogram(ideogram_bed_filename, chromosome, highlighted_region=None, color_map = giemsa_color_map, alignment="vertical", ax = None):
	ax = ax if ax is not None else plt.gca()

	ideogram_bed_file = open(ideogram_bed_filename, "r")

	# load fragments into list
	fragment_list_parm = []
	fragment_list_qarm = []
	for line in ideogram_bed_file:
		if(line[0] == "#"):
			continue

		split_line = line.rstrip().split("\t")

		if(not(split_line[0] == chromosome)):
			continue

		chrom = split_line[0]
		start = int(split_line[1])
		end = int(split_line[2])
		fragment_id = split_line[3]
		giemsa_stain = split_line[4]

		if(fragment_id[0] == "p"):
			fragment_list_parm += [[chrom, start, end, fragment_id, giemsa_stain]]
		elif(fragment_id[0] == "q"):
			fragment_list_qarm += [[chrom, start, end, fragment_id, giemsa_stain]]
	
	ideogram_bed_file.close()

	p_start = fragment_list_parm[0][1] if not len(fragment_list_parm) == 0 else None
	p_end = fragment_list_parm[-1][2] if not len(fragment_list_parm) == 0 else None
	q_start = fragment_list_qarm[0][1] if not len(fragment_list_qarm) == 0 else None
	q_end = fragment_list_qarm[-1][2] if not len(fragment_list_qarm) == 0 else None

	cap_ext = 2000000

	clip_path_patch_parm = None
	if(not p_start is None):
		# Create Clip path
		if(alignment == "vertical"):
			vertices = [[0, -1*cap_ext], [0, -1*p_end+cap_ext], [.5, bezierMidPoint(-1.*p_end, -1*p_end+cap_ext, -1*p_end+cap_ext, 0.5)], [1, -1*p_end+cap_ext], [1, -1.*cap_ext], [.5, bezierMidPoint(0, -1.*cap_ext, -1.*cap_ext, .5)], [0, -1*cap_ext]]
		else:
			vertices = [[1*cap_ext, 0], [1*p_end-cap_ext, 0], [bezierMidPoint(1.*p_end, 1*p_end-cap_ext, 1*p_end-cap_ext, 0.5), .5], [1*p_end-cap_ext, 1], [1.*cap_ext, 1], [bezierMidPoint(0, 1.*cap_ext, 1.*cap_ext, .5), .5], [1*cap_ext, 0]]

		codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3]
		clip_path = Path(vertices, codes=codes)
		clip_path_patch_parm = patches.PathPatch(clip_path, facecolor="none", linewidth=.5, capstyle="round")
		ax.add_patch(clip_path_patch_parm)

		rect = None
		if(alignment == "vertical"):
			rect = Rectangle((0, -1*p_end), 1, p_end, facecolor="none")
		else:
			rect = Rectangle((1*p_end, 0), p_end, 1, facecolor="none")
		rect.set_clip_path(clip_path_patch_parm)
		ax.add_patch(rect)

		# Plot fragments
		plotFragments(fragment_list_parm, clip_path_patch_parm, color_map, alignment=alignment, ax=ax)
	clip_path_path_qarm = None
	if(not q_start is None):
		if(alignment == "vertical"):
			vertices = [[0, -1*q_start-cap_ext], [0, -1*q_end+cap_ext], [.5, bezierMidPoint(-1*q_end, -1*q_end+cap_ext, -1*q_end+cap_ext, .5)], [1, -1*q_end+cap_ext], [1, -1*q_start-cap_ext], [.5, bezierMidPoint(-1*q_start, -1*q_start-cap_ext, -1*q_start-cap_ext, .5)], [0, -1*q_start-cap_ext]]
		else:
			vertices = [[1*q_start+cap_ext, 0], [1*q_end-cap_ext, 0], [bezierMidPoint(1*q_end, 1*q_end-cap_ext, 1*q_end-cap_ext, .5), .5], [1*q_end-cap_ext, 1], [1*q_start+cap_ext, 1], [bezierMidPoint(1*q_start, 1*q_start+cap_ext, 1*q_start+cap_ext, .5), .5], [1*q_start+cap_ext, 0]]

		codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3]
		clip_path = Path(vertices, codes=codes)
		clip_path_patch_qarm = patches.PathPatch(clip_path, facecolor="none", linewidth=.5, capstyle="round")
		ax.add_patch(clip_path_patch_qarm)

		rect = None
		if(alignment == "vertical"):
			rect = Rectangle((1, -1*q_end), 1, (q_end-q_start), facecolor="none")
		else:
			rect = Rectangle((1*q_end, 1), (q_end-q_start), 1, facecolor="none")
		rect.set_clip_path(clip_path_patch_qarm)
		ax.add_patch(rect)

		# Plot fragments
		plotFragments(fragment_list_qarm, clip_path_patch_qarm, color_map, alignment=alignment, ax=ax)

	ax.axis('off')

	if(alignment == "vertical"):
		plt.xlim([-.2, 1.2])
		plt.ylim([-1*q_end, 0])
	else:
		plt.ylim([-.2, 1.2])
		plt.xlim([0, 1*q_end])

	if(not highlighted_region is None):
		if (alignment == "horizontal"):
			rect = Rectangle((highlighted_region[1], 0), width=highlighted_region[2]-highlighted_region[1], height=1, capstyle="butt", linewidth=1, edgecolor="r", facecolor="none")
			ax.add_patch(rect)
