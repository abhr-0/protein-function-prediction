{
  description = "Protein function prediction: human Swiss-Prot (Track A) + CAFA6 (Track B)";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-26.05";
    flake-parts.url = "github:hercules-ci/flake-parts";
  };

  outputs =
    inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      perSystem =
        { pkgs, lib, ... }:
        let
          # nixpkgs only pulls in zlib/bzip2 for the Clang build. The pipeline
          # feeds MMseqs2 an uncompressed FASTA (run_clustering normalizes it
          # first), so this override is for parity with the Clang build and
          # ad-hoc use on compressed inputs, not a pipeline requirement.
          mmseqs2 = pkgs.mmseqs2.overrideAttrs (old: {
            buildInputs = (old.buildInputs or [ ]) ++ [
              pkgs.zlib
              pkgs.bzip2
            ];
            # Without a git checkout, MMseqs2 records "GITDIR-NOTFOUND" as its
            # version, which ends up in the split provenance JSON.
            cmakeFlags = (old.cmakeFlags or [ ]) ++ [
              (lib.cmakeFeature "VERSION_OVERRIDE" old.version)
            ];
          });
        in
        {
          devShells.default = pkgs.mkShell {
            packages = [
              pkgs.python312
              pkgs.uv
              mmseqs2
            ];

            env = lib.optionalAttrs pkgs.stdenv.isLinux {
              # Python wheels load native shared objects with dlopen(3); make the
              # loader aware of the manylinux libraries without RPATH lookups.
              LD_LIBRARY_PATH = lib.makeLibraryPath pkgs.pythonManylinuxPackages.manylinux1;
              # Use the Nix-provided interpreter instead of downloading one.
              UV_PYTHON_DOWNLOADS = "never";
            };

            shellHook = ''
              unset PYTHONPATH
            '';
          };
        };
    };
}
