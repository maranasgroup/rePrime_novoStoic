package nonEssential;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;


import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;

import com.opencsv.CSVReader;

public class organizeEmolecules {

	public static void main(String[] args) {		
		try {
			Connection conn = Connect.getConnection();
			Statement st = conn.createStatement();
			String database = "emolecules";
			String query = "select distinct smiles from "+database+".prices where SMILES is not null AND NEWSMILES is null and QUOTED_PRICE is not null";
			ResultSet rs = st.executeQuery(query);
			String smilesOriginal;
			String newSmiles;
			Statement updtStmnt = conn.createStatement();
			while (rs.next()) {
				smilesOriginal = rs.getString("smiles");
				newSmiles = getChemaxoneSmiles(smilesOriginal);
				if (null != newSmiles) {
					try {
						query = "UPDATE " + database + ".prices set NEWSMILES = '" + newSmiles + "' WHERE SMILES = '" + smilesOriginal + "'";
						//updtStmnt.executeUpdate(query);
						updtStmnt.addBatch(query);
					}
					catch (SQLException ex) {
						System.err.println(ex.getMessage());
					}
				}
			}
			updtStmnt.executeBatch();
			conn.close();
			st.close();
		}
		catch (SQLException ex) {
			System.err.println(ex.getMessage());
		}
	}

	private static String getChemaxoneSmiles(String smilesOriginal) {
		try {
			Molecule mol = MolImporter.importMol(smilesOriginal);
			Hydrogenize.convertExplicitHToImplicit(mol);		
			Cleaner.clean(mol, 2);  
			 
			String returnVal = MolExporter.exportToFormat(mol, "smiles:u0");
			String[] split = returnVal.split("[.]");
			int length = 0;
			for (int i = 0; i < split.length; i++) {
				if (split[i].length() > length) {
					returnVal = split[i]; 
				}
			}
			return returnVal;
		}
		catch (MolFormatException e) {
			System.out.println("error in "+smilesOriginal + " "+ e.getMessage());			
		}
		catch (IOException e) {
			System.out.println("error in "+smilesOriginal + " "+ e.getMessage());			
		}
		return null;
	}

}
