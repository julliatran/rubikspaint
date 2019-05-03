Project's name: Rubiks Paint

Description: Rubiks Paint allows the user to paint an imported .jpg image 
using only rubiks cube. The program will show the user how the cube should be
rotated with a 3D cube on screen, rotating the solution. The user can hit spacebar
to loop through the moves or select the speed to slow down to easily follow the 
solutions.

In addition to the main feature, there is also an open cube mode that allows the
user to just interact with the cube using the keyboard keys. Instructions are
displayed on screen.When the user import the image, a new file of the image 
is created using only rubiks cube color, effectively creating a pixel art style 
image from their picture. 

Libraries: To run this, the user need to install panda3d library, numpy, PIL and scipy.misc
    To download panda3d:
        1. Use this link https://www.panda3d.org/download.php?sdk
        
        2. Click 1.9.4 (the Latest stable) package and install. 
    
    To download PIL:     
        1. My PIL library is 1.1.7, install onto the python2.7 of panda3d, not 
        to be confused with the newer python version you might have. This PIL 
        library must be installed into python2.7 inorder for panda3d to import 
        properly so that imageToRubiks.py can run. 
        Here is a reference link: 
        https://stackoverflow.com/questions/20060096/installing-pil-with-pip
        module-manager didn't work for me so I had to use brew to install and
        terminal on mac to installed this. In python2.7, PIL is Pillow.

        2. For Mac, I used this command in terminal: 
            python2.7 -m pip install Pillow
        followed by: 
            sudo apt-get install libjpeg-dev zlib1g-dev 
        because it doesn't have that libjpeg library. 

        3. Make sure you also have pip installed because python2.7 doesn't have pip.

        To install pip on mac, you can do the command:
                    sudo easy_install pip 
        Here is a reference link: 
        https://stackoverflow.com/questions/50805485/cant-install-pip-on-python-2-7-only-python-3

    To download numpy:
        1. I installed using Homebrew and used the command:
            python2.7 pip install numpy
        Here is a reference link: 
        https://docs.scipy.org/doc/numpy-1.10.1/user/install.html#mac-os-x

    To install scipy:
        1. I installed using Homebrew with the command: 
            brew install gfortran
        2. Then: pip install scipy

        Make sure you have PIL installed first before you install scipy.
        Here is a reference link:
        https://penandpants.com/2012/02/24/install-python/

How to run: 1. Download the zip file, keeping everything in one folder. 
            2. Import needed libraries.
            3. Open main.py to run. 

Shortcut commads: No shortcut. To get to different screens, the user can just 
press the buttons because my screen display depends on user inputs. 
