# Docker_permission_denied

If you face ‘permission denied’ issues, you might need to add your user
to the Docker group. Follow these steps:

**1. Create the Docker Group (if it doesn’t already exist):**

``` sh
sudo groupadd docker
```

This command creates a Docker group if it does not already exist on your
system.

**2. Add Your User to the Docker Group:**

``` sh
sudo usermod -aG docker $USER
```

This command adds your user to the Docker group, which allows you to run
Docker commands without sudo.

**3. Log Out and Log Back In:**

This step is necessary to apply the group membership changes. After
logging back in, you should be able to run Docker commands without sudo.

Alternatively, you can use the following command to apply the group
changes immediately without logging out and back in:

``` sh
newgrp docker
```

**4. Verify Docker Permissions:**

To check if the permission issue is resolved, run a simple Docker
command:

``` sh
docker run hello-world
```

If the command runs successfully, you have configured Docker correctly.

**5. Ensure Docker Service is Running:**

If you still encounter issues, ensure that the Docker service is
running:

``` sh
sudo systemctl start docker
sudo systemctl enable docker
```

Following these steps will help you resolve permission issues and enable
you to pull and run the Docker container without using `sudo`.

## Session Info

    #> R version 4.5.2 (2025-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 22.04.5 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=de_AT.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: Europe/Vienna
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices datasets  utils     methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] desc_1.4.3          digest_0.6.39       R6_2.6.1           
    #>  [4] fastmap_1.2.0       xfun_0.56           cachem_1.1.0       
    #>  [7] knitr_1.51          htmltools_0.5.9     rmarkdown_2.30     
    #> [10] lifecycle_1.0.5     cli_3.6.5           sass_0.4.10        
    #> [13] pkgdown_2.2.0       textshaping_1.0.4   jquerylib_0.1.4    
    #> [16] renv_1.1.7          systemfonts_1.3.1   compiler_4.5.2     
    #> [19] rstudioapi_0.18.0   tools_4.5.2         ragg_1.5.0         
    #> [22] bslib_0.10.0        evaluate_1.0.5      yaml_2.3.12        
    #> [25] otel_0.2.0          BiocManager_1.30.27 jsonlite_2.0.0     
    #> [28] htmlwidgets_1.6.4   rlang_1.1.7         fs_1.6.6
