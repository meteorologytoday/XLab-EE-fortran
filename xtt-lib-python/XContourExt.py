import os;
import numpy as np;

def manualCLabelIfNotExists(ax, CS, file, folder='_clabel_info', rotation=True, **kwgs):
	file= "%s/%s.npy" % (folder, file,);

	try:
		os.mkdir(folder);
	except OSError:
		print('Directory \"%s\" already exists.'%(folder,));
		
	try:
		cl_pos = np.load(file);
		text_list = ax.clabel(CS, manual=cl_pos, **kwgs);
	except IOError as e:
		print("Clabel file not exists, use manual mode...");
		text_list = ax.clabel(CS, manual=True, **kwgs);
		pos = [_text.get_position() for _text in text_list];
		np.save(file, pos);
	
	if rotation == False:
		for key, value in enumerate(text_list):
			value.set_rotation(0);