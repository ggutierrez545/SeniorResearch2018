# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:22:48 2017

@author: gjgut
"""

import SavedGraphs
from moviepy.editor import *
import imageio
imageio.plugins.ffmpeg.download()

InfectGif = ImageSequenceClip('SavedGraphs', fps = 1,load_images = True)
InfectGif.write_gif('Inf_sem2.gif')