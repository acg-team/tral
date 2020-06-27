# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://vagrantcloud.com/search.
  config.vm.box = "ubuntu/bionic64"

  # Increase CPUs and memory for the box
  # https://stackoverflow.com/a/37335639/81658
  config.vm.provider "virtualbox" do |v|
    host = RbConfig::CONFIG['host_os']
    # Give VM 3/4 system memory & access to all cpu cores on the host
    if host =~ /darwin/
      cpus = `sysctl -n hw.ncpu`.to_i
    elsif host =~ /linux/
      cpus = `nproc`.to_i
    else # Windows folks
      cpus = `wmic cpu get NumberOfCores`.split("\n")[2].to_i
    end

    puts "Provisioning VM with #{cpus} CPU"
    v.customize ["modifyvm", :id, "--cpus", cpus]
  end

  # Need extra disk size for casper compilation
  # Requires disksize plugin:
  #   vagrant plugin install vagrant-disksize
  config.disksize.size = '15GB'

  # Log in as root
  #config.ssh.username = 'root'
  #config.ssh.password = 'vagrant'
  #config.ssh.insert_key = 'true'

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  # Enable provisioning with a shell script. Additional provisioners such as
  # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
  # documentation for more information about their specific syntax and use.
  # config.vm.provision "shell", inline: <<-SHELL
  #   apt-get update
  #   apt-get install -y apache2
  # SHELL
  config.vm.provision "shell", inline: <<-SHELL
    set -e

    # Install dependencies
    apt-get update
    apt-get -y dist-upgrade
    apt-get -y install python3 python3-pip unzip openjdk-8-jre-headless cmake
    pip3 install virtualenv

    # reduce wget output during provisioning
    echo 'verbose = off' >> ~/.wgetrc

    # Assistive tech cause java problems
    rm -f /usr/lib/jvm/java-8-openjdk-amd64/jre/lib/accessibility.properties

    cd /vagrant/easy_setup

    # accept all licenses
    if grep -q ACCEPT_ALL configTRAL_path.cfg; then
      sed -i 's/ACCEPT_ALL=no/ACCEPT_ALL=yes/i' configTRAL_path.cfg
    else
      echo ACCEPT_ALL=yes >> configTRAL_path.cfg
    fi


    # Install TRAL software
    ./setupTRAL.sh setup

    # Config file
    cat <<END > ~/.tral/config.ini
###########################################
### Configuration file for TRAL Vagrant ###
###########################################

sequence_type = AA

[sequence]
    [[repeat_detection]]
        # AA includes all detectors used by default on protein sequence data.
        AA = HHrepID, T-REKS, TRUST, XSTREAM
        # DNA includes all detectors used by default on protein sequence data.
        DNA = Phobos, TRED, T-REKS, TRF, XSTREAM
    [[repeat_detector_path]]
        # If the executable is in the system path, supply its name. Otherwise, supply the full path to the executable. Details are explained in TRAL's online docs.
        PHOBOS = phobos
        HHrepID = hhrepid_64
        HHrepID_dummyhmm = ~/.tral/data/hhrepid/dummyHMM.hmm
        T-REKS = T-REKS
        TRED = tred
        TRF = trf
        TRUST = TRUST
        TRUST_substitutionmatrix = ~/.tral/tral_external_software/TRUST_Align/Align/BLOSUM50
        XSTREAM = XSTREAM

[hmm]
    hmmbuild = hmmbuild
    l_effective_max = 50

[filter]
    [[basic]]
        tag = basic_filter
        [[[dict]]]
            [[[[pvalue]]]]
                func_name = pvalue
                score = phylo_gap01
                threshold = 0.1
            [[[[n_effective]]]]
                func_name = attribute
                attribute = n_effective
                type = min
                threshold = 1.9


[repeat]
    scoreslist = phylo_gap01, # score (the comma in the end is needed for TRAL)
    calc_score = False # is the score calculated?
    calc_pvalue = False # is the pvalue calculated?
    precision = 10
    ginsi = ginsi # integrated in MAFFT
    Castor = Castor
    [[castor_parameter]]
        rate_distribution = constant # either constant or gamma
    alfsim = alfsim

[repeat_list]
    # Columns to include in repeat list TSV output
    # Allowed values:
    # - begin: position of the tandem repeats within the sequence,
    # - pvalue: statistical significance of the tandem repeats
    # - divergence: divergence of the tandem repeat units
    # - l_effective: length of the tandem repeat units
    # - n_effective: number of tandem repeat units
    # - msa_original: multiple sequence alignment
    # - score: score corresponding to the value of 'model'
    # - repeat_region_length: total length of repeat region
    output_characteristics = begin, msa_original, l_effective, n_effective, repeat_region_length, divergence, pvalue

    # model for scoring repeats. Supported: entropy, parsimony, pSim, phylo, phylo_gap01, phylo_gap001
    model = phylo_gap01

[repeat_score]
    evolutionary_model = lg
    [[indel]]
        indel_rate_per_site = 0.01
        ignore_gaps = True
        gaps = row_wise
        zipf = 1.821
    [[optimisation]]
        start_min = 0.5
        start_max = 1.5
        n_iteration = 14
    [[K80]]
        kappa = 2.59
    [[TN93]]
        alpha_1 = 0.3
        alpha_2 = 0.4
        beta = 0.7
    [[score_calibration]]
        scoreslist=phylo_gap01, # score (the comma at the end is needed)
        save_calibration = False
        precision = 10

[AA]
    standard_chars = A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    all_chars = A, B, C, D, E, F, G, H, I, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    [[ambiguous_chars]]
        B = D,N
        O = K,
        U = C,
        Z = E,Q
        X = A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
[DNA]
    standard_chars = A, C, G, T
    all_chars = A, C, G, T, N, X
    [[ambiguous_chars]]
        N = A, C, G, T
        X = A, C, G, T

END


    # All external software
    ./install_ext_software.sh

    cd /vagrant
    # dev requirements are optional but useful for tests and docs
    pip3 install -r requirements_dev.txt

    echo
    echo "THIS MACHINE CONTAINS PROPRIETARY SOFTWARE."
    echo "Please check the licenses before using (e.g. no commercial use permitted)"
  SHELL

end
