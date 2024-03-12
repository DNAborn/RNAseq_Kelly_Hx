update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     \
                    libblas.so-x86_64-linux-gnu      /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so 50
                    

sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     \
                    libblas.so-x86_64-linux-gnu      /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so 50
sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3   \
                    libblas.so.3-x86_64-linux-gnu    /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so 50
sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   \
                    liblapack.so-x86_64-linux-gnu    /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so 50
sudo update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 \
                    liblapack.so.3-x86_64-linux-gnu  /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so 50
                    # 
                    
echo "/opt/intel/lib/intel64"     >  /etc/ld.so.conf.d/mkl.conf
echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf
ldconfig

sudo sh -c "echo '/opt/intel/oneapi/2024.0/lib'     >  /etc/ld.so.conf.d/mkl.conf"
sudo sh -c "echo '/opt/intel/oneapi/mkl/latest/lib' >> /etc/ld.so.conf.d/mkl.conf"
sudo ldconfig

find /opt/intel -name "libmkl_sycl.so"

sudo rm /opt/intel/oneapi/2024.0/lib/libmkl_sycl.so
sudo rm /opt/intel/oneapi/mkl/2024.0/lib/libmkl_sycl.so

sudo ldconfig

sudo sh -c "echo 'MKL_THREADING_LAYER=GNU' >> /etc/environment"

sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu
sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu
