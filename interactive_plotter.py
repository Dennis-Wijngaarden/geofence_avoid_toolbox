#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 14:33:19 2019

@author: dennis
"""
import numpy as np
import matplotlib.pyplot as plt
import linear_tools as lin_tools

x_geo_0 = 400.
x_geo_1 = 500.
y_geo_0 = -500.
y_geo_1 = 500.

a_geo, b_geo, c_geo = lin_tools.line_eq_from_points((x_geo_0, y_geo_0), (x_geo_1, y_geo_1))
d_int_x = 100.0
d_int_y = 1.0
d_own_x = -500.0
d_own_y = -0.0
v_int_x = -1.0
v_int_y = -0.0
v_res_x = np.arange(-20., 20., 0.001)

v_res_y_1 = (-a_geo*d_int_y*v_res_x + 2*a_geo*d_own_x*v_int_y + a_geo*d_own_y*v_res_x + b_geo*d_int_x*v_int_x - b_geo*d_int_x*v_res_x + b_geo*d_int_y*v_int_y - b_geo*d_own_x*v_int_x + b_geo*d_own_x*v_res_x + b_geo*d_own_y*v_int_y + 2*c_geo*v_int_y - np.sqrt(4*a_geo**2*d_int_x*d_own_x*v_int_x*v_res_x - 4*a_geo**2*d_int_x*d_own_x*v_res_x**2 + a_geo**2*d_int_y**2*v_res_x**2 - 2*a_geo**2*d_int_y*d_own_y*v_res_x**2 - 4*a_geo**2*d_own_x**2*v_int_x**2 + 4*a_geo**2*d_own_x**2*v_int_x*v_res_x + a_geo**2*d_own_y**2*v_res_x**2 + 2*a_geo*b_geo*d_int_x*d_int_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_x*d_int_y*v_res_x**2 + 4*a_geo*b_geo*d_int_x*d_own_x*v_int_x*v_int_y - 4*a_geo*b_geo*d_int_x*d_own_x*v_int_y*v_res_x + 2*a_geo*b_geo*d_int_x*d_own_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_x*d_own_y*v_res_x**2 + 2*a_geo*b_geo*d_int_y**2*v_int_y*v_res_x - 4*a_geo*b_geo*d_int_y*d_own_x*v_int_x**2 + 6*a_geo*b_geo*d_int_y*d_own_x*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_y*d_own_x*v_res_x**2 - 4*a_geo*b_geo*d_int_y*d_own_y*v_int_y*v_res_x - 4*a_geo*b_geo*d_own_x**2*v_int_x*v_int_y + 4*a_geo*b_geo*d_own_x**2*v_int_y*v_res_x - 4*a_geo*b_geo*d_own_x*d_own_y*v_int_x**2 + 6*a_geo*b_geo*d_own_x*d_own_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_own_x*d_own_y*v_res_x**2 + 2*a_geo*b_geo*d_own_y**2*v_int_y*v_res_x + 4*a_geo*c_geo*d_int_x*v_int_x*v_res_x - 4*a_geo*c_geo*d_int_x*v_res_x**2 - 8*a_geo*c_geo*d_own_x*v_int_x**2 + 12*a_geo*c_geo*d_own_x*v_int_x*v_res_x - 4*a_geo*c_geo*d_own_x*v_res_x**2 + b_geo**2*d_int_x**2*v_int_x**2 - 2*b_geo**2*d_int_x**2*v_int_x*v_res_x + b_geo**2*d_int_x**2*v_res_x**2 + 2*b_geo**2*d_int_x*d_int_y*v_int_x*v_int_y - 2*b_geo**2*d_int_x*d_int_y*v_int_y*v_res_x - 2*b_geo**2*d_int_x*d_own_x*v_int_x**2 + 4*b_geo**2*d_int_x*d_own_x*v_int_x*v_res_x - 2*b_geo**2*d_int_x*d_own_x*v_res_x**2 + 2*b_geo**2*d_int_x*d_own_y*v_int_x*v_int_y - 2*b_geo**2*d_int_x*d_own_y*v_int_y*v_res_x + b_geo**2*d_int_y**2*v_int_y**2 - 2*b_geo**2*d_int_y*d_own_x*v_int_x*v_int_y + 2*b_geo**2*d_int_y*d_own_x*v_int_y*v_res_x - 4*b_geo**2*d_int_y*d_own_y*v_int_x**2 + 8*b_geo**2*d_int_y*d_own_y*v_int_x*v_res_x - 2*b_geo**2*d_int_y*d_own_y*v_int_y**2 - 4*b_geo**2*d_int_y*d_own_y*v_res_x**2 + b_geo**2*d_own_x**2*v_int_x**2 - 2*b_geo**2*d_own_x**2*v_int_x*v_res_x + b_geo**2*d_own_x**2*v_res_x**2 - 2*b_geo**2*d_own_x*d_own_y*v_int_x*v_int_y + 2*b_geo**2*d_own_x*d_own_y*v_int_y*v_res_x + b_geo**2*d_own_y**2*v_int_y**2 + 4*b_geo*c_geo*d_int_x*v_int_x*v_int_y - 4*b_geo*c_geo*d_int_x*v_int_y*v_res_x - 4*b_geo*c_geo*d_int_y*v_int_x**2 + 8*b_geo*c_geo*d_int_y*v_int_x*v_res_x - 4*b_geo*c_geo*d_int_y*v_res_x**2 - 4*b_geo*c_geo*d_own_x*v_int_x*v_int_y + 4*b_geo*c_geo*d_own_x*v_int_y*v_res_x - 4*b_geo*c_geo*d_own_y*v_int_x**2 + 8*b_geo*c_geo*d_own_y*v_int_x*v_res_x - 4*b_geo*c_geo*d_own_y*v_res_x**2 - 4*c_geo**2*v_int_x**2 + 8*c_geo**2*v_int_x*v_res_x - 4*c_geo**2*v_res_x**2))/(2*(a_geo*d_own_x + b_geo*d_int_y + c_geo))
v_res_y_2 = (-a_geo*d_int_y*v_res_x + 2*a_geo*d_own_x*v_int_y + a_geo*d_own_y*v_res_x + b_geo*d_int_x*v_int_x - b_geo*d_int_x*v_res_x + b_geo*d_int_y*v_int_y - b_geo*d_own_x*v_int_x + b_geo*d_own_x*v_res_x + b_geo*d_own_y*v_int_y + 2*c_geo*v_int_y + np.sqrt(4*a_geo**2*d_int_x*d_own_x*v_int_x*v_res_x - 4*a_geo**2*d_int_x*d_own_x*v_res_x**2 + a_geo**2*d_int_y**2*v_res_x**2 - 2*a_geo**2*d_int_y*d_own_y*v_res_x**2 - 4*a_geo**2*d_own_x**2*v_int_x**2 + 4*a_geo**2*d_own_x**2*v_int_x*v_res_x + a_geo**2*d_own_y**2*v_res_x**2 + 2*a_geo*b_geo*d_int_x*d_int_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_x*d_int_y*v_res_x**2 + 4*a_geo*b_geo*d_int_x*d_own_x*v_int_x*v_int_y - 4*a_geo*b_geo*d_int_x*d_own_x*v_int_y*v_res_x + 2*a_geo*b_geo*d_int_x*d_own_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_x*d_own_y*v_res_x**2 + 2*a_geo*b_geo*d_int_y**2*v_int_y*v_res_x - 4*a_geo*b_geo*d_int_y*d_own_x*v_int_x**2 + 6*a_geo*b_geo*d_int_y*d_own_x*v_int_x*v_res_x - 2*a_geo*b_geo*d_int_y*d_own_x*v_res_x**2 - 4*a_geo*b_geo*d_int_y*d_own_y*v_int_y*v_res_x - 4*a_geo*b_geo*d_own_x**2*v_int_x*v_int_y + 4*a_geo*b_geo*d_own_x**2*v_int_y*v_res_x - 4*a_geo*b_geo*d_own_x*d_own_y*v_int_x**2 + 6*a_geo*b_geo*d_own_x*d_own_y*v_int_x*v_res_x - 2*a_geo*b_geo*d_own_x*d_own_y*v_res_x**2 + 2*a_geo*b_geo*d_own_y**2*v_int_y*v_res_x + 4*a_geo*c_geo*d_int_x*v_int_x*v_res_x - 4*a_geo*c_geo*d_int_x*v_res_x**2 - 8*a_geo*c_geo*d_own_x*v_int_x**2 + 12*a_geo*c_geo*d_own_x*v_int_x*v_res_x - 4*a_geo*c_geo*d_own_x*v_res_x**2 + b_geo**2*d_int_x**2*v_int_x**2 - 2*b_geo**2*d_int_x**2*v_int_x*v_res_x + b_geo**2*d_int_x**2*v_res_x**2 + 2*b_geo**2*d_int_x*d_int_y*v_int_x*v_int_y - 2*b_geo**2*d_int_x*d_int_y*v_int_y*v_res_x - 2*b_geo**2*d_int_x*d_own_x*v_int_x**2 + 4*b_geo**2*d_int_x*d_own_x*v_int_x*v_res_x - 2*b_geo**2*d_int_x*d_own_x*v_res_x**2 + 2*b_geo**2*d_int_x*d_own_y*v_int_x*v_int_y - 2*b_geo**2*d_int_x*d_own_y*v_int_y*v_res_x + b_geo**2*d_int_y**2*v_int_y**2 - 2*b_geo**2*d_int_y*d_own_x*v_int_x*v_int_y + 2*b_geo**2*d_int_y*d_own_x*v_int_y*v_res_x - 4*b_geo**2*d_int_y*d_own_y*v_int_x**2 + 8*b_geo**2*d_int_y*d_own_y*v_int_x*v_res_x - 2*b_geo**2*d_int_y*d_own_y*v_int_y**2 - 4*b_geo**2*d_int_y*d_own_y*v_res_x**2 + b_geo**2*d_own_x**2*v_int_x**2 - 2*b_geo**2*d_own_x**2*v_int_x*v_res_x + b_geo**2*d_own_x**2*v_res_x**2 - 2*b_geo**2*d_own_x*d_own_y*v_int_x*v_int_y + 2*b_geo**2*d_own_x*d_own_y*v_int_y*v_res_x + b_geo**2*d_own_y**2*v_int_y**2 + 4*b_geo*c_geo*d_int_x*v_int_x*v_int_y - 4*b_geo*c_geo*d_int_x*v_int_y*v_res_x - 4*b_geo*c_geo*d_int_y*v_int_x**2 + 8*b_geo*c_geo*d_int_y*v_int_x*v_res_x - 4*b_geo*c_geo*d_int_y*v_res_x**2 - 4*b_geo*c_geo*d_own_x*v_int_x*v_int_y + 4*b_geo*c_geo*d_own_x*v_int_y*v_res_x - 4*b_geo*c_geo*d_own_y*v_int_x**2 + 8*b_geo*c_geo*d_own_y*v_int_x*v_res_x - 4*b_geo*c_geo*d_own_y*v_res_x**2 - 4*c_geo**2*v_int_x**2 + 8*c_geo**2*v_int_x*v_res_x - 4*c_geo**2*v_res_x**2))/(2*(a_geo*d_own_x + b_geo*d_int_y + c_geo))

nul_x_1 = (2*a_geo**2*d_int_x*d_own_x*v_int_x + 2*a_geo**2*d_own_x**2*v_int_x + a_geo*b_geo*d_int_x*d_int_y*v_int_x - 2*a_geo*b_geo*d_int_x*d_own_x*v_int_y + a_geo*b_geo*d_int_x*d_own_y*v_int_x + a_geo*b_geo*d_int_y**2*v_int_y + 3*a_geo*b_geo*d_int_y*d_own_x*v_int_x - 2*a_geo*b_geo*d_int_y*d_own_y*v_int_y + 2*a_geo*b_geo*d_own_x**2*v_int_y + 3*a_geo*b_geo*d_own_x*d_own_y*v_int_x + a_geo*b_geo*d_own_y**2*v_int_y + 2*a_geo*c_geo*d_int_x*v_int_x + 6*a_geo*c_geo*d_own_x*v_int_x - b_geo**2*d_int_x**2*v_int_x - b_geo**2*d_int_x*d_int_y*v_int_y + 2*b_geo**2*d_int_x*d_own_x*v_int_x - b_geo**2*d_int_x*d_own_y*v_int_y + b_geo**2*d_int_y*d_own_x*v_int_y + 4*b_geo**2*d_int_y*d_own_y*v_int_x - b_geo**2*d_own_x**2*v_int_x + b_geo**2*d_own_x*d_own_y*v_int_y - 2*b_geo*c_geo*d_int_x*v_int_y + 4*b_geo*c_geo*d_int_y*v_int_x + 2*b_geo*c_geo*d_own_x*v_int_y + 4*b_geo*c_geo*d_own_y*v_int_x + 4*c_geo**2*v_int_x - 2*np.sqrt((a_geo*d_own_x + b_geo*d_int_y + c_geo)*(a_geo*d_own_x + b_geo*d_own_y + c_geo)*(d_int_x**2 - 2*d_int_x*d_own_x + d_int_y**2 - 2*d_int_y*d_own_y + d_own_x**2 + d_own_y**2))*(a_geo*v_int_x + b_geo*v_int_y))/(4*a_geo**2*d_int_x*d_own_x - a_geo**2*d_int_y**2 + 2*a_geo**2*d_int_y*d_own_y - a_geo**2*d_own_y**2 + 2*a_geo*b_geo*d_int_x*d_int_y + 2*a_geo*b_geo*d_int_x*d_own_y + 2*a_geo*b_geo*d_int_y*d_own_x + 2*a_geo*b_geo*d_own_x*d_own_y + 4*a_geo*c_geo*d_int_x + 4*a_geo*c_geo*d_own_x - b_geo**2*d_int_x**2 + 2*b_geo**2*d_int_x*d_own_x + 4*b_geo**2*d_int_y*d_own_y - b_geo**2*d_own_x**2 + 4*b_geo*c_geo*d_int_y + 4*b_geo*c_geo*d_own_y + 4*c_geo**2)
nul_x_2 = (2*a_geo**2*d_int_x*d_own_x*v_int_x + 2*a_geo**2*d_own_x**2*v_int_x + a_geo*b_geo*d_int_x*d_int_y*v_int_x - 2*a_geo*b_geo*d_int_x*d_own_x*v_int_y + a_geo*b_geo*d_int_x*d_own_y*v_int_x + a_geo*b_geo*d_int_y**2*v_int_y + 3*a_geo*b_geo*d_int_y*d_own_x*v_int_x - 2*a_geo*b_geo*d_int_y*d_own_y*v_int_y + 2*a_geo*b_geo*d_own_x**2*v_int_y + 3*a_geo*b_geo*d_own_x*d_own_y*v_int_x + a_geo*b_geo*d_own_y**2*v_int_y + 2*a_geo*c_geo*d_int_x*v_int_x + 6*a_geo*c_geo*d_own_x*v_int_x - b_geo**2*d_int_x**2*v_int_x - b_geo**2*d_int_x*d_int_y*v_int_y + 2*b_geo**2*d_int_x*d_own_x*v_int_x - b_geo**2*d_int_x*d_own_y*v_int_y + b_geo**2*d_int_y*d_own_x*v_int_y + 4*b_geo**2*d_int_y*d_own_y*v_int_x - b_geo**2*d_own_x**2*v_int_x + b_geo**2*d_own_x*d_own_y*v_int_y - 2*b_geo*c_geo*d_int_x*v_int_y + 4*b_geo*c_geo*d_int_y*v_int_x + 2*b_geo*c_geo*d_own_x*v_int_y + 4*b_geo*c_geo*d_own_y*v_int_x + 4*c_geo**2*v_int_x + 2*np.sqrt((a_geo*d_own_x + b_geo*d_int_y + c_geo)*(a_geo*d_own_x + b_geo*d_own_y + c_geo)*(d_int_x**2 - 2*d_int_x*d_own_x + d_int_y**2 - 2*d_int_y*d_own_y + d_own_x**2 + d_own_y**2))*(a_geo*v_int_x + b_geo*v_int_y))/(4*a_geo**2*d_int_x*d_own_x - a_geo**2*d_int_y**2 + 2*a_geo**2*d_int_y*d_own_y - a_geo**2*d_own_y**2 + 2*a_geo*b_geo*d_int_x*d_int_y + 2*a_geo*b_geo*d_int_x*d_own_y + 2*a_geo*b_geo*d_int_y*d_own_x + 2*a_geo*b_geo*d_own_x*d_own_y + 4*a_geo*c_geo*d_int_x + 4*a_geo*c_geo*d_own_x - b_geo**2*d_int_x**2 + 2*b_geo**2*d_int_x*d_own_x + 4*b_geo**2*d_int_y*d_own_y - b_geo**2*d_own_x**2 + 4*b_geo*c_geo*d_int_y + 4*b_geo*c_geo*d_own_y + 4*c_geo**2)

#######
# GUI #
#######

class GeoPlotter(object):
    def __init__(self, client):
        self.client = client
        self.press = None
        self.background = None
        self.selected_arg = None
        self.selected_element = None
        
        self.fig = plt.figure()
        self.ax = plt.subplot()
        self.fig.canvas.set_window_title('top view')
        plt.grid()
        self.geofence, = self.ax.plot(self.client.geofence[:,0], self.client.geofence[:,1], label='geofence')
        self.p_own, = self.ax.plot(self.client.p_own[0], self.client.p_own[1], 'go', label='own ac')
        self.p_int, = self.ax.plot(self.client.p_int[0], self.client.p_int[1], 'ro', label='int ac')
    
    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
    
    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if not event.inaxes: return
        # check to which point is clicked closest
        dist_geo = self.calc_pixel_dist(self.geofence, event)
        dist_own = self.calc_pixel_dist(self.p_own, event)
        dist_int = self.calc_pixel_dist(self.p_int, event)
        
        dist_array = np.concatenate((dist_geo, dist_own, dist_int))
        self.selected_arg = np.argmin(dist_array)
        
        # Draw eevrything but the selected element and store the pixel buffer
        if (self.selected_arg < len(dist_geo)):
            self.selected_element = self.geofence
            x0 = self.geofence.get_data()[0][self.selected_arg]
            y0 = self.geofence.get_data()[1][self.selected_arg]
        elif (self.selected_arg == len(dist_geo)):
            self.selected_element = self.p_own
            x0 = self.p_own.get_data()[0][0]
            y0 = self.p_own.get_data()[1][0]
        else:
            self.selected_element = self.p_int
            x0 = self.p_int.get_data()[0][0]
            y0 = self.p_int.get_data()[1][0]
        
        self.press = x0, y0, event.xdata, event.ydata

        canvas = self.selected_element.figure.canvas
        axes = self.selected_element.axes
        self.selected_element.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.selected_element.axes.bbox)
        
        # now redraw selected element
        axes.draw_artist(self.selected_element)
        
        # Blit the redrawn area
        canvas.blit(axes.bbox)
    
    def on_motion(self, event):
        'On motion we will move the element'
        if self.selected_element == None: return
        if not event.inaxes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        
        if self.selected_element == self.geofence:
            x_array = self.geofence.get_data()[0]
            y_array = self.geofence.get_data()[1]
            
            x_array[self.selected_arg] = x0 + dx
            y_array[self.selected_arg] = y0 + dy
            
            self.selected_element.set_xdata(x_array)
            self.selected_element.set_ydata(y_array)
            self.client.geofence = np.array([x_array, y_array]).T
        else:
            x = x0 + dx
            y = y0 + dy
            self.selected_element.set_xdata(np.array([x]))
            self.selected_element.set_ydata(np.array([y]))
            
            if self.selected_element == self.p_own:
                self.client.p_own = np.array([x, y])
            else:
                self.client.p_int = np.array([x, y])
        
        canvas = self.selected_element.figure.canvas
        axes = self.selected_element.axes
        # restore the background region
        canvas.restore_region(self.background)
        
        # redraw just the selected element
        axes.draw_artist(self.selected_element)
        
        # blit just the selected element
        canvas.blit(axes.bbox)
        
    def on_release(self, event):
        'on release we reset the press data'
        if self.selected_element == None: return
        
        self.press = None
        
        # turn of the animation property, draw and reset background
        self.selected_element.set_animated(False)
        
        self.selected_element.figure.canvas.draw()
        
        self.selected_element = None
        self.background = None
        
    def calc_pixel_dist(self, points, event):
        xref, yref = event.x, event.y
        
        x, y = points.get_data()
        xy_pixels = self.ax.transData.transform(np.vstack([x,y]).T)
        xpix, ypix = xy_pixels.T
        
        dist = np.sqrt((xref-xpix)**2 + (yref-ypix)**2)
        return np.array(dist)
    
class SpeedPlotter(object):
    def __init__(self, client):
        self.client = client
        self.press = None
        self.background = None
        self.lock = False
        
        self.fig = plt.figure()
        self.ax = plt.subplot()
        self.fig.canvas.set_window_title('intruder speed')
        self.ax.set_xlim(-30., 30.)
        self.ax.set_ylim(-30., 30.)
        plt.grid()
        self.v_int, = self.ax.plot([0., self.client.v_int[0]], [0., self.client.v_int[1]], label='intruder speed')
        
    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        
    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if not event.inaxes: return
        self.lock = True
        
        # Draw eevrything but the selected element and store the pixel buffer
        x0 = self.v_int.get_data()[0][1]
        y0 = self.v_int.get_data()[1][1]
        
        self.press = x0, y0, event.xdata, event.ydata
        
        canvas = self.v_int.figure.canvas
        axes = self.v_int.axes
        self.v_int.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.v_int.axes.bbox)
        
        # now redraw selected element
        axes.draw_artist(self.v_int)
        
        # Blit the redrawn area
        canvas.blit(axes.bbox)
    
    def on_motion(self, event):
        'On motion we will move the element'
        if not self.lock: return
        if not event.inaxes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        
        x_array = np.array([0., x0 + dx])
        y_array = np.array([0., y0 + dy])
        
        self.v_int.set_xdata(x_array)
        self.v_int.set_ydata(y_array)
        self.client.v_int = np.array([x_array[1], y_array[1]])

        
        canvas = self.v_int.figure.canvas
        axes = self.v_int.axes
        # restore the background region
        canvas.restore_region(self.background)
        
        # redraw just the selected element
        axes.draw_artist(self.v_int)
        
        # blit just the selected element
        canvas.blit(axes.bbox)
    
    def on_release(self, event):
        'on release we reset the press data'
        if not self.lock: return
        
        self.press = None
        
        # turn of the animation property, draw and reset background
        self.v_int.set_animated(False)
        
        self.v_int.figure.canvas.draw()
        
        self.background = None
        self.lock = False
          
class Client(object):
    def __init__(self, geofence, p_own, p_int, v_int):
        # Triogoniometric parameters
        self.geofence = geofence
        self.p_own = p_own
        self.p_int = p_int
        self.v_int = v_int
        
        # Interactive plotters
        self.geo_plotter = GeoPlotter(self)
        self.geo_plotter.connect()
        
        self.speed_plotter = SpeedPlotter(self)
        self.speed_plotter.connect()

plt.close("all")
geofence = np.array([[x_geo_0, y_geo_0], [x_geo_1, y_geo_1]])
p_own = np.array([d_own_x, d_own_y])
p_int = np.array([d_int_x, d_int_y])
v_int = np.array([v_int_x, v_int_y])

client = Client(geofence, p_own, p_int, v_int)

############################
# Plotter for the geofence #
############################
