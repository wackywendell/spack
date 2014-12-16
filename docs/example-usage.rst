
Example Usage of spack
======================

This is an ipython notebook in ``spack/docs/example-usage.ipynb``,
demonstrating how to use ``spack``.

First, we import our modules:

.. code:: python

    import spack
    from math import pi
Now we put in some data.

The data below is from a simple packing I made using ``pyparm.packmin``.
Normally you'd load your own data from a file.

.. code:: python

    L = 2.0066668050219723
    diameters = [ 0.96,  0.97,  0.98,  0.99,  1.  ,  1.01,  1.02,  1.03,  1.04]
    locs = [[ 1.40776762,  1.26647724,  0.73389219],
            [ 0.58704249,  2.11399   ,  1.52956579],
            [ 1.75917911,  0.54290089,  1.27577478],
            [ 2.13750384,  0.87508242,  0.21938647],
            [ 1.07283961,  0.87692084,  1.9060841 ],
            [ 0.09550267,  1.94404465,  0.56463369],
            [ 1.07636871,  2.1942971 ,  0.63752152],
            [ 0.49922725,  1.20002224,  1.13360082],
            [-0.27724757,  1.62152603,  1.67262247]]
Now we make the packing:

.. code:: python

    pack = spack.Packing(locs, diameters, L=L)
What does it look like?

.. code:: python

    size = 400
    sc = pack.scene(rot=pi/4, camera_dist=2, cmap='autumn', bgcolor=[1,1,1])
    sc.render('ipython', width=size, height=size, antialiasing=0.001)



.. image:: example-usage_files/example-usage_9_0.png



Let's make a movie!

.. code:: python

    from moviepy.editor import VideoClip
    import moviepy.editor as mpy
    
    duration = 10
    def make_frame(t):
        return (
                    pack.scene(rot=(t/duration + .125)*2.*pi, 
                               camera_dist=2, cmap='autumn', bgcolor=[1,1,1])
                    .render(width=size, height=size, antialiasing=0.001)
                )
    vc = VideoClip(make_frame, duration=duration)
Write the movie to file. This takes about 10 minutes on my machine to
render all 240 frames, but it does give you some pretty spiffy output.

.. code:: python

    vc.write_gif("example-packing.gif",fps=24)

.. parsed-literal::

    
    [MoviePy] >>>> Building file example-packing.gif
    [MoviePy] Generating GIF frames...
    [MoviePy] Optimizing the GIF with ImageMagick...
    [MoviePy] >>>> File example-packing.gif is ready !

And display the movie, in an ipython notebook.

.. code:: python

    from IPython.display import Image
    Image(url="example-packing.gif")



.. raw:: html

    <img src="example-packing.gif"/>


