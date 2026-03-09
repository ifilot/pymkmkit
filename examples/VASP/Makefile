# Directory containing example data
DATA_DIR := ../../tests/data

# Archives
ZIP_NI := $(DATA_DIR)/OUTCAR_Ni311_C.zip
ZIP_RU := $(DATA_DIR)/OUTCAR_Ru1121_C.zip
ZIP_TS := $(DATA_DIR)/OUTCAR_Ru1121_C_CH_TS.zip

# Extracted OUTCAR files
OUTCAR_NI := $(DATA_DIR)/OUTCAR_Ni311_C
OUTCAR_RU := $(DATA_DIR)/OUTCAR_Ru1121_C
OUTCAR_TS := $(DATA_DIR)/OUTCAR_Ru1121_C_CH_TS

# YAML outputs
YAML_NI := ni311_c.yaml
YAML_RU := ru1121_c.yaml
YAML_TS := ru1121_c_ch_ts.yaml

.PHONY: all clean show

# --------------------------------------------------
# Build all examples
# --------------------------------------------------

all: $(YAML_NI) $(YAML_RU) $(YAML_TS)
	@echo ""
	@echo "Generated example YAML files:"
	@echo "  $(YAML_NI)"
	@echo "  $(YAML_RU)"
	@echo "  $(YAML_TS)"
	@echo ""

# --------------------------------------------------
# Extraction
# --------------------------------------------------

$(OUTCAR_NI): $(ZIP_NI)
	@echo "Extracting Ni OUTCAR..."
	unzip -o $< -d $(DATA_DIR)

$(OUTCAR_RU): $(ZIP_RU)
	@echo "Extracting Ru OUTCAR (paired modes)..."
	unzip -o $< -d $(DATA_DIR)

$(OUTCAR_TS): $(ZIP_TS)
	@echo "Extracting Ru TS OUTCAR..."
	unzip -o $< -d $(DATA_DIR)

# --------------------------------------------------
# YAML generation
# --------------------------------------------------

$(YAML_NI): $(OUTCAR_NI)
	@echo "Generating YAML (Ni system)..."
	pymkmkit freq2yaml $< -o $@

$(YAML_RU): $(OUTCAR_RU)
	@echo "Generating YAML (Ru double-adsorbate)..."
	pymkmkit freq2yaml $< -o $@ --average-pairs

$(YAML_TS): $(OUTCAR_TS)
	@echo "Generating YAML (Ru transition state)..."
	pymkmkit freq2yaml $< -o $@ --average-pairs

# --------------------------------------------------
# Show results
# --------------------------------------------------

show: $(YAML_NI) $(YAML_RU) $(YAML_TS)
	@echo ""
	@echo "===== $(YAML_NI) ====="
	@cat $(YAML_NI)
	@echo ""
	@echo "===== $(YAML_RU) ====="
	@cat $(YAML_RU)
	@echo ""
	@echo "===== $(YAML_TS) ====="
	@cat $(YAML_TS)

# --------------------------------------------------
# Clean generated files
# --------------------------------------------------

clean:
	rm -f \
	$(OUTCAR_NI) \
	$(OUTCAR_RU) \
	$(OUTCAR_TS) \
	$(YAML_NI) \
	$(YAML_RU) \
	$(YAML_TS)
