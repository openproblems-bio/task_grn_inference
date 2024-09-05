# Use the base image
FROM continuumio/anaconda3:2024.02-1

# Install Python packages
RUN pip install anndata

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz
RUN tar -xvzf julia-1.9.3-linux-x86_64.tar.gz
RUN mv julia-1.9.3 /opt/
RUN ln -s /opt/julia-1.9.3/bin/julia /usr/local/bin/julia

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.resolve(); Pkg.add(PackageSpec(name="InformationMeasures", version="0.3.0")); Pkg.add(PackageSpec(name="PyPlot", version="2.8.0")); Pkg.add(PackageSpec(name="LightGraphs", version="1.2.0")); Pkg.add(PackageSpec(name="GraphPlot", version="0.3.1")); Pkg.add(PackageSpec(name="NetworkInference", version="0.1.0"))'
