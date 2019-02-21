#!/bin/bash

# visualize relative abundances of taxonomic groups with ampvis2

# get mothur .wang.taxonomy files
cp ../*.wang.taxonomy .
for i in *.both.*.taxonomy; do mv $i ${i/both/1000bp}; done

# convert taxonomy format to something ampvis2 can read
for i in *.wang.taxonomy; do

    # take taxonomy paths, introduce k__ at start, and replace ; at the end with a dummy s__
    # mothur classify.seqs() didn't report a 7th field, so creating a dummy one here
    cut -f 2 $i | sed -r "s/\([0-9]+\)//g" | sed -r -e "s/^/k__/" -e "s/\;$/\|s__/" > $i.ampvis2

    # convert taxonomy path format
    sed -i -r -e "s/\;/|p__/1" -e "s/\;/|c__/1" -e "s/\;/|o__/1" -e "s/\;/|f__/1" -e "s/\;/|g__/1" $i.ampvis2

    # count occurrence of each taxonomic path
    sort $i.ampvis2 | uniq -c | sed -r "s/^\s+([0-9]+)\s/\1\t/" > $i.ampvis2.count

    # introduce header
    sample=$(echo $i | cut -f 2 -d .)
    length=$(echo $i | cut -f 5 -d . )
    sed -i "1s/^/${sample}_${length}\tKingdom|Phylum|Class|Order|Family|Genus|Species\n/" $i.ampvis2.count
done

# merge files
## sequential joining of tables
join -a 1 -a 2 -e '0' -1 2 -2 2 -o '1.1,2.1,0' -t $'\t' \
    rnammer.P19.prok.16S.1000bp.nr_v128.wang.taxonomy.ampvis2.count rnammer.P19.prok.16S.250bp.nr_v128.wang.taxonomy.ampvis2.count > dummy1
join -a 1 -a 2 -e '0' -1 2 -2 3 -o '1.1,2.1,2.2,0' -t $'\t' \
    rnammer.PM3.prok.16S.250bp.nr_v128.wang.taxonomy.ampvis2.count dummy1 > dummy2
join -a 1 -a 2 -e '0' -1 2 -2 4 -o '1.1,2.1,2.2,2.3,0' -t $'\t' \
    rnammer.PM3.prok.16S.1000bp.nr_v128.wang.taxonomy.ampvis2.count dummy2 > dummy3
join -a 1 -a 2 -e '0' -1 2 -2 5 -o '1.1,2.1,2.2,2.3,2.4,0' -t $'\t' \
    rnammer.SALA.prok.16S.250bp.nr_v128.wang.taxonomy.ampvis2.count dummy3 > dummy4
join -a 1 -a 2 -e '0' -1 2 -2 6 -o '1.1,2.1,2.2,2.3,2.4,2.5,0' -t $'\t' \
    rnammer.SALA.prok.16S.1000bp.nr_v128.wang.taxonomy.ampvis2.count dummy4 > dummy5
join -a 1 -a 2 -e '0' -1 2 -2 7 -o '1.1,2.1,2.2,2.3,2.4,2.5,2.6,0' -t $'\t' \
    rnammer.TNS08.prok.16S.250bp.nr_v128.wang.taxonomy.ampvis2.count dummy5 > dummy6
join -a 1 -a 2 -e '0' -1 2 -2 8 -o '1.1,2.1,2.2,2.3,2.4,2.5,2.6,2.7,0' -t $'\t' \
    rnammer.TNS08.prok.16S.1000bp.nr_v128.wang.taxonomy.ampvis2.count dummy6 > allSamples.ampvis2.count
rm dummy*
# introduce dummy OTU column
paste <( seq 324 | sed -r "s/^/OTU_/" | sed "1s/^/OTU\n/" ) allSamples.ampvis2.count | sed -r "s/\|/\t/g" > allSamples.ampvis2.count.final

# in R do

archaea.class<-c("Aigarchaeota; Terrestrial_Hot_Spring_Gp(THSCG)","Aigarchaeota; Aigarchaeota_Incertae_Sedis","Aigarchaeota; Aigarchaeota_unclassified","Bathyarchaeota; Bathyarchaeota_cl","Crenarchaeota; Thermoprotei","Euryarchaeota; Methanomicrobia","Euryarchaeota; Thermoplasmata","Euryarchaeota; Euryarchaeota_unclassified","Hadesarchaea; Hadesarchaea_cl","Thaumarchaeota; AK59","Thaumarchaeota; Group_C3","Thaumarchaeota; Marine_Benthic_Group_A","Thaumarchaeota; Marine_Group_I","Thaumarchaeota; OPPD003","Thaumarchaeota; South_African_Gold_Mine_Gp_1(SAGMCG-1)","Thaumarchaeota; Thaumarchaeota_unclassified","Thaumarchaeota; YS18As63","Woesearchaeota_(DHVEG-6); Woesearchaeota_(DHVEG-6)_cl","Marine_Hydrothermal_Vent_Group(MHVG); Marine_Hydrothermal_Vent_Group(MHVG)_cl","pMC2A209; pMC2A209_cl","Candidate_division_YNPFFA; Candidate_division_YNPFFA_cl","WSA2; 19c-33","AK8; AK8_cl","Archaea_unclassified; Archaea_unclassified")
bacteria.class<-c("Acidobacteria; Acidobacteria","Acidobacteria; Blastocatellia","Acidobacteria; Subgroup_17","Acidobacteria; Subgroup_19","Acidobacteria; Subgroup_22","Acidobacteria; Subgroup_6","Acidobacteria; Acidobacteria_unclassified","Actinobacteria; Acidimicrobiia","Actinobacteria; Actinobacteria","Actinobacteria; MB-A2-108","Actinobacteria; OPB41","Actinobacteria; Thermoleophilia","Actinobacteria; Actinobacteria_unclassified","Aerophobetes; Aerophobetes_cl","Aminicenantes; Aminicenantes_cl","Aquificae; Aquificae","Armatimonadetes; uncultured","Atribacteria; Atribacteria_cl","Bacteroidetes; Bacteroidetes_BD2-2","Bacteroidetes; Bacteroidetes_vadinHA17","Bacteroidetes; Bacteroidia","Bacteroidetes; Cytophagia","Bacteroidetes; Flavobacteriia","Bacteroidetes; SB-5","Bacteroidetes; Sphingobacteriia","Bacteroidetes; Bacteroidetes_unclassified","Chlamydiae; LD1-PA32","Chlorobi; Chlorobia","Chloroflexi; Anaerolineae","Chloroflexi; Ardenticatenia","Chloroflexi; Caldilineae","Chloroflexi; Chloroflexi_cl","Chloroflexi; JG30-KF-CM66","Chloroflexi; KD4-96","Chloroflexi; MSB-5E12","Chloroflexi; S085","Chloroflexi; SBR2076","Chloroflexi; SJA-68","Chloroflexi; Thermoflexia","Chloroflexi; Thermomicrobia","Chloroflexi; TK10","Chloroflexi; uncultured","Chloroflexi; Chloroflexi_unclassified","Cyanobacteria; Chloroplast","Cyanobacteria; Cyanobacteria","Cyanobacteria; Cyanobacteria_unclassified","Deferribacteres; Deferribacteres","Deferribacteres; Deferribacteres_Incertae_Sedis","Dictyoglomi; Dictyoglomia","Elusimicrobia; Elusimicrobia","Fervidibacteria; Fervidibacteria_cl","Firmicutes; Bacilli","Firmicutes; Clostridia","Firmicutes; Negativicutes","Gemmatimonadetes; BD2-11_terrestrial_group","Gemmatimonadetes; Gemmatimonadetes","Gemmatimonadetes; MD2902-B12","Gemmatimonadetes; Gemmatimonadetes_unclassified","Gracilibacteria; Gracilibacteria_cl","Hydrogenedentes; Hydrogenedentes_cl","Hydrogenedentes; Unknown_Class","Ignavibacteriae; Ignavibacteria","Latescibacteria; Latescibacteria_cl","Latescibacteria; Latescibacteria_Incertae_Sedis","Latescibacteria; Latescibacteria_unclassified","Marinimicrobia_(SAR406_clade); Marinimicrobia_(SAR406_clade)_cl","Nitrospirae; Nitrospira","Omnitrophica; Omnitrophica_cl","Omnitrophica; Omnitrophica_Incertae_Sedis","Omnitrophica; Omnitrophica_unclassified","Planctomycetes; MD2896-B258","Planctomycetes; OM190","Planctomycetes; Pla4_lineage","Planctomycetes; Planctomycetacia","Planctomycetes; vadinHA49","Planctomycetes; Planctomycetes_unclassified","Proteobacteria; Alphaproteobacteria","Proteobacteria; Betaproteobacteria","Proteobacteria; Gammaproteobacteria","Proteobacteria; Deltaproteobacteria","Proteobacteria; Epsilonproteobacteria","Proteobacteria; MACA-EFT26","Proteobacteria; SPOTSOCT00m83","Proteobacteria; Proteobacteria_unclassified","Spirochaetae; Spirochaetes","Tectomicrobia; Tectomicrobia_cl","Tenericutes; Mollicutes","Thermodesulfobacteria; Thermodesulfobacteria","Thermotogae; Thermotogae","Verrucomicrobia; OPB35_soil_group","Verrucomicrobia; Opitutae","Verrucomicrobia; Verrucomicrobia_Incertae_Sedis","Verrucomicrobia; Verrucomicrobia_unclassified","RBG-1_(Zixibacteria); RBG-1_(Zixibacteria)_cl","SBR1093; SBR1093_cl","TM6_(Dependentiae); TM6_(Dependentiae)_cl","WS1; WS1_cl","WS2; WS2_cl","GAL15; GAL15_cl","Bacteria_unclassified; Bacteria_unclassified")
unclassified<-c("unknown_unclassified; unknown_unclassified")
order.class<-c(archaea.class,bacteria.class,unclassified)

library("ampvis2")
myotutable<-read.table("allSamples.ampvis2.count.final", check.names=F, header=T, row.names=1)
mymetadata<-read.table("metadata.csv", header=T, colClasses="character")
d <- amp_load(otutable = myotutable, metadata=mymetadata)
amp_heatmap(d, tax_show=50, plot_na=F, tax_aggregate = "Class", tax_add = "Phylum", group_by = "Length", facet_by = "Name", order_y_by=rev(order.class), plot_values_size=3) + theme(axis.text.x = element_text(angle=45, size=9, vjust=1), axis.text.y = element_text(size=9))
ggsave(file="compareRelAbund.pdf", width = 210, height = 260, units = "mm")

