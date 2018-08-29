### mafft: Alignment of tandem repeat units
LINK_MAFFT=https://mafft.cbrc.jp/alignment/software/mafft_7.397-1_amd64.deb
sudo wget $LINK_MAFFT -P $TRAL_TOOLS
sudo dpkg -i $TRAL_TOOLS/mafft_7.397-1_amd64.deb
sudo rm -rf $TRAL_TOOLS/mafft_7.397-1_amd64.deb
