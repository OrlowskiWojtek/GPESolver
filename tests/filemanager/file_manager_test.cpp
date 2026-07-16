#include <gtest/gtest.h>
#include "filemanager/file_manager.hpp"
#include "nlohmann/json.hpp"

// Fixtures for common test setup
class FileManagerTest : public ::testing::Test {
protected:
    AbstractSimulationMediator* mediator;
    FileManager* fileManager;
    nlohmann::json validJson;

    void SetUp() override {
        mediator = nullptr;
        fileManager = new FileManager(mediator);
        
        // Valid base configuration
        validJson = {
            {"n_atoms", 1000},
            {"m", 4.0},
            {"dx", 0.1},
            {"dy", 0.1},
            {"dz", 0.1},
            {"nx", 64},
            {"ny", 64},
            {"nz", 64},
            {"omega_x", 1.0e9},
            {"omega_y", 1.0e9},
            {"omega_z", 1.0e9},
            {"calc_strategy", "IT"},
            {"init_strategy", "GAUSS"},
            {"pote_strategy", "HARMONIC"},
            {"iter_imag", 1000},
            {"iter_real", 5000},
            {"edd", 0.0},
            {"fftw_n_threads", 4}
        };
    }

    void TearDown() override {
        delete fileManager;
    }
};

// ============================================
// load_box tests
// ============================================

TEST_F(FileManagerTest, LoadBox_AllRequiredParams_Success) {
    // Given: validJson has all required box params (dx, dy, dz, nx, ny, nz)
    nlohmann::json j = validJson;

    // When/Then: should not throw
    EXPECT_NO_THROW(fileManager->load_box(j));
}

TEST_F(FileManagerTest, LoadBox_MissingDx_Throws) {
    nlohmann::json j = validJson;
    j.erase("dx");

    EXPECT_THROW(fileManager->load_box(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadBox_MissingNx_Throws) {
    nlohmann::json j = validJson;
    j.erase("nx");

    EXPECT_THROW(fileManager->load_box(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadBox_MissingAllRequired_Throws) {
    nlohmann::json j;
    // Empty json - all required params missing

    EXPECT_THROW(fileManager->load_box(j), std::runtime_error);
}

// ============================================
// load_simulation tests
// ============================================

TEST_F(FileManagerTest, LoadSimulation_AllRequiredParams_Success) {
    nlohmann::json j = validJson;

    EXPECT_NO_THROW(fileManager->load_simulation(j));
}

TEST_F(FileManagerTest, LoadSimulation_MissingCalcStrategy_Throws) {
    nlohmann::json j = validJson;
    j.erase("calc_strategy");

    EXPECT_THROW(fileManager->load_simulation(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadSimulation_MissingNIterImag_Throws) {
    nlohmann::json j = validJson;
    j.erase("iter_imag");

    EXPECT_THROW(fileManager->load_simulation(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadSimulation_MissingNAtoms_Throws) {
    nlohmann::json j = validJson;
    j.erase("n_atoms");

    EXPECT_THROW(fileManager->load_simulation(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadSimulation_MissingMass_Throws) {
    nlohmann::json j = validJson;
    j.erase("m");

    EXPECT_THROW(fileManager->load_simulation(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadSimulation_OptionalFftwThreads_UsesDefault) {
    nlohmann::json j = validJson;
    j.erase("fftw_n_threads");

    // Should not throw, and use default value (4)
    EXPECT_NO_THROW(fileManager->load_simulation(j));
}

// ============================================
// load_initialization tests
// ============================================

TEST_F(FileManagerTest, LoadInitialization_SingleGauss_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "GAUSS";
    j.erase("initial_maximas");

    EXPECT_NO_THROW(fileManager->load_initialization(j));
}

TEST_F(FileManagerTest, LoadInitialization_MultipleGauss_MissingInitialMaximas_Throws) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "MULTIPLE_GAUSS";
    j.erase("initial_maximas");
    // Missing "initial_maximas"

    EXPECT_THROW(fileManager->load_initialization(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadInitialization_MultipleGauss_WithInitialMaximas_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "MULTIPLE_GAUSS";
    j["initial_maximas"] = 5;

    EXPECT_NO_THROW(fileManager->load_initialization(j));
}

TEST_F(FileManagerTest, LoadInitialization_SetupGauss_MissingBecDroplets_Throws) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "setup_gauss";
    // Missing bec_droplets_x, y, z

    EXPECT_THROW(fileManager->load_initialization(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadInitialization_SetupGauss_WithAllBecDroplets_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "SETUP_GAUSS";
    j["bec_droplets_x"] = 1.0;
    j["bec_droplets_y"] = 1.0;
    j["bec_droplets_z"] = 1.0;

    EXPECT_NO_THROW(fileManager->load_initialization(j));
}

TEST_F(FileManagerTest, LoadInitialization_FromBinaryFile_MissingLoadFilename_Throws) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "BINARY_FILE";
    // Missing "load_filename"

    EXPECT_THROW(fileManager->load_initialization(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadInitialization_WrongInitStrategy) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "WRONG_NAME";

    EXPECT_THROW(fileManager->load_initialization(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadInitialization_FromBinaryFile_WithLoadFilename_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "BINARY_FILE";
    j["load_filename"] = "test.gpe.bin";

    EXPECT_NO_THROW(fileManager->load_initialization(j));
}

TEST_F(FileManagerTest, LoadInitialization_FromTextFile_WithLoadFilename_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "TEXT_FILE";
    j["load_filename"] = "test.gpe.dat";

    EXPECT_NO_THROW(fileManager->load_initialization(j));
}

// ============================================
// load_potential tests
// ============================================

TEST_F(FileManagerTest, LoadPotential_Harmonic_Success) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "HARMONIC";

    EXPECT_NO_THROW(fileManager->load_potential(j));
}

TEST_F(FileManagerTest, LoadPotential_Mexican_MissingDd_Throws) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "MEXICAN";
    j.erase("dd");  // Remove dd - should throw

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadPotential_Mexican_WithDd_Success) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "MEXICAN";
    j["dd"] = 0.5;

    EXPECT_NO_THROW(fileManager->load_potential(j));
}

TEST_F(FileManagerTest, LoadPotential_MexicanFree_MissingDd_Throws) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "MEXICAN_FREE";
    j.erase("dd");

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadPotential_MexicanAsymetric_MissingDd_Throws) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "MEXICAN_ASYMETRIC";
    j.erase("dd");

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadPotential_MissingOmega_Throws) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "HARMONIC";
    j.erase("omega_x");

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadPotential_MissingPoteStrategy_Throws) {
    nlohmann::json j = validJson;
    j.erase("pote_strategy");

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

TEST_F(FileManagerTest, LoadPotential_WrongPoteStrategy_Throws) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "IAMNOTREAL";

    EXPECT_THROW(fileManager->load_potential(j), std::runtime_error);
}

// ============================================
// load_all_v1 tests (integration test)
// ============================================

TEST_F(FileManagerTest, LoadAllV1_CompleteConfig_Success) {
    nlohmann::json j = validJson;

    EXPECT_NO_THROW(fileManager->load_all_v1(j));
}

TEST_F(FileManagerTest, LoadAllV1_MultipleGaussScenario_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "MULTIPLE_GAUSS";
    j["initial_maximas"] = 3;

    EXPECT_NO_THROW(fileManager->load_all_v1(j));
}

TEST_F(FileManagerTest, LoadAllV1_SetupGaussScenario_Success) {
    nlohmann::json j = validJson;
    j["init_strategy"] = "SETUP_GAUSS";
    j["bec_droplets_x"] = 1.0;
    j["bec_droplets_y"] = 1.0;
    j["bec_droplets_z"] = 1.0;
    j["pote_strategy"] = "MEXICAN";
    j["dd"] = 0.5;

    EXPECT_NO_THROW(fileManager->load_all_v1(j));
}

TEST_F(FileManagerTest, LoadAllV1_MexicanPotentialRequiresDd_Fails) {
    nlohmann::json j = validJson;
    j["pote_strategy"] = "MEXICAN";
    // Missing dd - should fail because load_potential requires it for MEXICAN

    EXPECT_THROW(fileManager->load_all_v1(j), std::runtime_error);
}

// ============================================
// Error message tests
// ============================================

TEST_F(FileManagerTest, CheckRequired_ErrorMessageContainsKey) {
    nlohmann::json j;
    j["other_key"] = 42;

    try {
        fileManager->load_box(j);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("dx"), std::string::npos);
    }
}
