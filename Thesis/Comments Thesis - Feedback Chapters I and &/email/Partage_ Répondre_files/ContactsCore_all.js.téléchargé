if (AjxPackage.define("ContactsCore")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2006, 2007, 2008, 2009, 2010, 2012, 2013, 2014, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2006, 2007, 2008, 2009, 2010, 2012, 2013, 2014, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */
/*
 * Package: ContactsCore
 * 
 * Supports: Loading of contacts and address books
 * 
 * Loaded:
 * 	- When user contacts are loaded during startup
 * 	- If the <refresh> block has address books
 * 	- If a search for contacts returns results
 */

if (AjxPackage.define("zimbraMail.abook.model.ZmAddrBook")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */

/**
 * @overview
 * This file contains the address book class.
 */

/**
 * Creates an address book.
 * @constructor
 * @class
 * This class represents an address book.
 * 
 * @author Parag Shah
 *
 * @param	{Hash}	params		a hash of parameters
 * @param {int}	params.id			a numeric ID
 * @param {String}	params.name		[string]	the name
 * @param {ZmOrganizer}	params.parent		the parent organizer
 * @param {ZmTree}	params.tree		a tree model that contains this organizer
 * @param {int}	params.color		the color of this address book
 * @param {String}	params.owner	the owner of the address book (if shared)
 * @param {String}	params.zid 		the the share ID of a shared address book
 * @param {String}	params.rid		the the remote folder id of a shared address book
 * @param {String}	params.restUrl	the REST URL of this organizer
 * 
 * @extends		ZmFolder
 */
ZmAddrBook = function(params) {
	params.type = ZmOrganizer.ADDRBOOK;
	ZmFolder.call(this, params);
};

ZmAddrBook.prototype = new ZmFolder;
ZmAddrBook.prototype.constructor = ZmAddrBook;

ZmAddrBook.prototype.isZmAddrBook = true;
ZmAddrBook.prototype.toString = function() { return "ZmAddrBook"; };


// Consts

ZmAddrBook.ID_ADDRESSBOOK = ZmOrganizer.ID_ADDRBOOK; 							// XXX: may not be necessary

// Public methods

ZmAddrBook.prototype.getIcon =
function() {
	if (this.nId == ZmFolder.ID_ROOT)			{ return null; }
	if (this.nId == ZmFolder.ID_TRASH)			{ return "Trash"; }
	if (this.link || this.isRemote())			{ return "SharedContactsFolder"; }
	if (this.nId == ZmFolder.ID_AUTO_ADDED)		{ return "EmailedContacts"; }
	return "ContactsFolder";
};

/**
 * Checks if the address book supports public access.
 * 
 * @return	{Boolean}		always returns <code>true</code>
 */
ZmAddrBook.prototype.supportsPublicAccess =
function() {
	// AddrBook's can be accessed outside of ZCS (i.e. REST)
	return true;
};

/**
 * @private
 */
ZmAddrBook.prototype.mayContain = function(what) {

	if (!what) {
		return true;
	}

	// Distribution Lists is a system-generated folder
	if (this.id == ZmOrganizer.ID_DLS) {
		return false;
	}

	if (what.isZmAddrBook) {
		return ZmFolder.prototype.mayContain.apply(this, arguments);
	}

	// An item or an array of items is being moved
	var items = AjxUtil.toArray(what);
	for (var i = 0; i < items.length; i++) {
		var item = items[i];
		if (item.type !== ZmItem.CONTACT && item.type !== ZmItem.GROUP) {
			// only contacts are valid for addr books.
			return false;
		}
	}

	return ZmFolder.prototype.mayContain.apply(this, arguments);
};
}
if (AjxPackage.define("zimbraMail.abook.model.ZmContact")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */

/**
 * @overview
 * This file contains the contact class.
 */

if (!window.ZmContact) {
/**
 * Creates an empty contact.
 * @class
 * This class represents a contact (typically a person) with all its associated versions
 * of email address, home and work addresses, phone numbers, etc. Contacts can be filed/sorted
 * in different ways, with the default being Last, First. A contact is an item, so
 * it has tagging and flagging support, and belongs to a list.
 * <p>
 * Most of a contact's data is kept in attributes. These include name, phone, etc. Meta-data and
 * data common to items are not kept in attributes. These include flags, tags, folder, and
 * modified/created dates. Since the attribute data for contacts is loaded only once, a contact
 * gets its attribute values from that canonical list.
 * </p>
 *
 * @param {int}	id		the unique ID
 * @param {ZmContactList}	list		the list that contains this contact
 * @param {constant}	type		the item type
 * @param {object}	newDl		true if this is a new DL
 *
 * @extends		ZmItem
 */
ZmContact = function(id, list, type, newDl) {
	if (arguments.length == 0) { return; }

	type = type || ZmItem.CONTACT;
	ZmItem.call(this, type, id, list);

	this.attr = {};
	this.isGal = (this.list && this.list.isGal) || newDl;
	if (newDl) {
		this.folderId = ZmFolder.ID_DLS;
		this.dlInfo = {	isMember: false,
						isOwner: true,
						subscriptionPolicy: null,
						unsubscriptionPolicy: null,
						description: "",
						displayName: "",
						notes: "",
						hideInGal: false,
						mailPolicy: null,
						owners: [appCtxt.get(ZmSetting.USERNAME)]
		};

	}

	this.participants = new AjxVector(); // XXX: need to populate this guy (see ZmConv)
};

ZmContact.prototype = new ZmItem;
ZmContact.prototype.constructor = ZmContact;
ZmContact.prototype.isZmContact = true;

// fields
ZmContact.F_anniversary				= "anniversary";
ZmContact.F_assistantPhone			= "assistantPhone";
ZmContact.F_attachment				= "attachment";
ZmContact.F_birthday				= "birthday";
ZmContact.F_callbackPhone			= "callbackPhone";
ZmContact.F_carPhone				= "carPhone";
ZmContact.F_company					= "company";
ZmContact.F_companyPhone			= "companyPhone";
ZmContact.F_custom					= "custom";
ZmContact.F_description				= "description";
ZmContact.F_department				= "department";
ZmContact.F_dlist					= "dlist";				// Group fields
ZmContact.F_dlDisplayName			= "dldisplayname"; //DL
ZmContact.F_dlDesc					= "dldesc";  //DL
ZmContact.F_dlHideInGal				= "dlhideingal";  //DL
ZmContact.F_dlNotes					= "dlnotes";  //DL
ZmContact.F_dlSubscriptionPolicy	= "dlsubspolicy";  //DL
ZmContact.F_dlMailPolicy			= "dlmailpolicy";  //DL
ZmContact.F_dlMailPolicySpecificMailers	= "dlmailpolicyspecificmailers";  //DL
ZmContact.F_dlUnsubscriptionPolicy	= "dlunsubspolicy";  //DL
ZmContact.F_dlListOwners			= "dllistowners";  //DL
ZmContact.F_email					= "email";
ZmContact.F_email2					= "email2";
ZmContact.F_email3					= "email3";
ZmContact.F_email4					= "email4";
ZmContact.F_email5					= "email5";
ZmContact.F_email6					= "email6";
ZmContact.F_email7					= "email7";
ZmContact.F_email8					= "email8";
ZmContact.F_email9					= "email9";
ZmContact.F_email10					= "email10";
ZmContact.F_email11					= "email11";
ZmContact.F_email12					= "email12";
ZmContact.F_email13					= "email13";
ZmContact.F_email14					= "email14";
ZmContact.F_email15					= "email15";
ZmContact.F_email16					= "email16";
ZmContact.F_fileAs					= "fileAs";
ZmContact.F_firstName				= "firstName";
ZmContact.F_folderId				= "folderId";
ZmContact.F_groups                  = "groups";         //group members
ZmContact.F_homeCity				= "homeCity";
ZmContact.F_homeCountry				= "homeCountry";
ZmContact.F_homeFax					= "homeFax";
ZmContact.F_homePhone				= "homePhone";
ZmContact.F_homePhone2				= "homePhone2";
ZmContact.F_homePostalCode			= "homePostalCode";
ZmContact.F_homeState				= "homeState";
ZmContact.F_homeStreet				= "homeStreet";
ZmContact.F_homeURL					= "homeURL";
ZmContact.F_image					= "image";				// contact photo
ZmContact.F_imAddress 				= "imAddress";			// IM addresses
ZmContact.F_imAddress1 				= "imAddress1";			// IM addresses
ZmContact.F_imAddress2 				= "imAddress2";
ZmContact.F_imAddress3				= "imAddress3";
ZmContact.F_jobTitle				= "jobTitle";
ZmContact.F_lastName				= "lastName";
ZmContact.F_maidenName				= "maidenName";
ZmContact.F_memberC                 = "memberC";
ZmContact.F_memberG                 = "memberG";
ZmContact.F_memberI                 = "memberI";
ZmContact.F_middleName				= "middleName";
ZmContact.F_mobilePhone				= "mobilePhone";
ZmContact.F_namePrefix				= "namePrefix";
ZmContact.F_nameSuffix				= "nameSuffix";
ZmContact.F_nickname				= "nickname";
ZmContact.F_notes					= "notes";
ZmContact.F_otherCity				= "otherCity";
ZmContact.F_otherCountry			= "otherCountry";
ZmContact.F_otherFax				= "otherFax";
ZmContact.F_otherPhone				= "otherPhone";
ZmContact.F_otherPostalCode			= "otherPostalCode";
ZmContact.F_otherState				= "otherState";
ZmContact.F_otherStreet				= "otherStreet";
ZmContact.F_otherURL				= "otherURL";
ZmContact.F_pager					= "pager";
ZmContact.F_phoneticFirstName       = "phoneticFirstName";
ZmContact.F_phoneticLastName        = "phoneticLastName";
ZmContact.F_phoneticCompany         = "phoneticCompany";
ZmContact.F_type					= "type";
ZmContact.F_workAltPhone			= "workAltPhone";
ZmContact.F_workCity				= "workCity";
ZmContact.F_workCountry				= "workCountry";
ZmContact.F_workEmail1				= "workEmail1";
ZmContact.F_workEmail2				= "workEmail2";
ZmContact.F_workEmail3				= "workEmail3";
ZmContact.F_workFax					= "workFax";
ZmContact.F_workMobile				= "workMobile";
ZmContact.F_workPhone				= "workPhone";
ZmContact.F_workPhone2				= "workPhone2";
ZmContact.F_workPostalCode			= "workPostalCode";
ZmContact.F_workState				= "workState";
ZmContact.F_workStreet				= "workStreet";
ZmContact.F_workURL					= "workURL";
ZmContact.F_imagepart               = "imagepart";          // New field for bug 73146 - Contacts call does not return the image information
ZmContact.F_zimletImage				= "zimletImage";
ZmContact.X_fileAs					= "fileAs";				// extra fields
ZmContact.X_firstLast				= "firstLast";
ZmContact.X_fullName				= "fullName";
ZmContact.X_vcardXProps             = "vcardXProps";
ZmContact.X_outlookUserField        = "outlookUserField";
ZmContact.MC_cardOwner				= "cardOwner";			// My card fields
ZmContact.MC_workCardMessage		= "workCardMessage";
ZmContact.MC_homeCardMessage		= "homeCardMessage";
ZmContact.MC_homePhotoURL			= "homePhotoURL";
ZmContact.MC_workPhotoURL			= "workPhotoURL";
ZmContact.GAL_MODIFY_TIMESTAMP		= "modifyTimeStamp";	// GAL fields
ZmContact.GAL_CREATE_TIMESTAMP		= "createTimeStamp";
ZmContact.GAL_ZIMBRA_ID				= "zimbraId";
ZmContact.GAL_OBJECT_CLASS			= "objectClass";
ZmContact.GAL_MAIL_FORWARD_ADDRESS	= "zimbraMailForwardingAddress";
ZmContact.GAL_CAL_RES_TYPE			= "zimbraCalResType";
ZmContact.GAL_CAL_RES_LOC_NAME		= "zimbraCalResLocationDisplayName";

// file as
(function() {
	var i = 1;
	ZmContact.FA_LAST_C_FIRST			= i++;
	ZmContact.FA_FIRST_LAST 			= i++;
	ZmContact.FA_COMPANY 				= i++;
	ZmContact.FA_LAST_C_FIRST_COMPANY	= i++;
	ZmContact.FA_FIRST_LAST_COMPANY		= i++;
	ZmContact.FA_COMPANY_LAST_C_FIRST	= i++;
	ZmContact.FA_COMPANY_FIRST_LAST		= i++;
	ZmContact.FA_CUSTOM					= i++;
})();

// Field information

ZmContact.ADDRESS_FIELDS = [
    // NOTE: sync with field order in ZmEditContactView's templates
	ZmContact.F_homeCity,
	ZmContact.F_homeCountry,
	ZmContact.F_homePostalCode,
	ZmContact.F_homeState,
	ZmContact.F_homeStreet,
	ZmContact.F_workCity,
	ZmContact.F_workCountry,
	ZmContact.F_workPostalCode,
	ZmContact.F_workState,
	ZmContact.F_workStreet,
    ZmContact.F_otherCity,
    ZmContact.F_otherCountry,
    ZmContact.F_otherPostalCode,
    ZmContact.F_otherState,
    ZmContact.F_otherStreet
];
ZmContact.EMAIL_FIELDS = [
	ZmContact.F_email,
	ZmContact.F_workEmail1,
	ZmContact.F_workEmail2,
	ZmContact.F_workEmail3
];
ZmContact.IM_FIELDS = [
	ZmContact.F_imAddress
];
ZmContact.OTHER_FIELDS = [
    // NOTE: sync with field order in ZmEditContactView's templates
	ZmContact.F_birthday,
    ZmContact.F_anniversary,
	ZmContact.F_custom
];
ZmContact.PHONE_FIELDS = [
    // NOTE: sync with field order in ZmEditContactView's templates
    ZmContact.F_mobilePhone,
    ZmContact.F_workPhone,
    ZmContact.F_workFax,
    ZmContact.F_companyPhone,
    ZmContact.F_homePhone,
    ZmContact.F_homeFax,
    ZmContact.F_pager,
    ZmContact.F_callbackPhone,
	ZmContact.F_assistantPhone,
	ZmContact.F_carPhone,
	ZmContact.F_otherPhone,
    ZmContact.F_otherFax,
	ZmContact.F_workAltPhone,
	ZmContact.F_workMobile
];
ZmContact.PRIMARY_FIELDS = [
    // NOTE: sync with field order in ZmEditContactView's templates
    ZmContact.F_image,
    ZmContact.F_namePrefix,
    ZmContact.F_firstName,
    ZmContact.F_phoneticFirstName,
    ZmContact.F_middleName,
	ZmContact.F_maidenName,
    ZmContact.F_lastName,
    ZmContact.F_phoneticLastName,
    ZmContact.F_nameSuffix,
    ZmContact.F_nickname,
    ZmContact.F_jobTitle,
    ZmContact.F_department,
	ZmContact.F_company,
    ZmContact.F_phoneticCompany,
	ZmContact.F_fileAs,
	ZmContact.F_folderId,
	ZmContact.F_notes
];
ZmContact.URL_FIELDS = [
    // NOTE: sync with field order in ZmEditContactView's templates
	ZmContact.F_homeURL,
	ZmContact.F_workURL,
	ZmContact.F_otherURL
];
ZmContact.GAL_FIELDS = [
	ZmContact.GAL_MODIFY_TIMESTAMP,
	ZmContact.GAL_CREATE_TIMESTAMP,
	ZmContact.GAL_ZIMBRA_ID,
	ZmContact.GAL_OBJECT_CLASS,
	ZmContact.GAL_MAIL_FORWARD_ADDRESS,
	ZmContact.GAL_CAL_RES_TYPE,
	ZmContact.GAL_CAL_RES_LOC_NAME,
	ZmContact.F_type
];
ZmContact.MYCARD_FIELDS = [
	ZmContact.MC_cardOwner,
	ZmContact.MC_homeCardMessage,
	ZmContact.MC_homePhotoURL,
	ZmContact.MC_workCardMessage,
	ZmContact.MC_workPhotoURL
];
ZmContact.X_FIELDS = [
	ZmContact.X_firstLast,
	ZmContact.X_fullName,
    ZmContact.X_vcardXProps
];


ZmContact.IGNORE_NORMALIZATION = [];

ZmContact.ADDR_PREFIXES = ["work","home","other"];
ZmContact.ADDR_SUFFIXES = ["Street","City","State","PostalCode","Country"];

ZmContact.updateFieldConstants = function() {

	for (var i = 0; i < ZmContact.ADDR_PREFIXES.length; i++) {
		for (var j = 0; j < ZmContact.ADDR_SUFFIXES.length; j++) {
			ZmContact.IGNORE_NORMALIZATION.push(ZmContact.ADDR_PREFIXES[i] + ZmContact.ADDR_SUFFIXES[j]);
		}
	}

ZmContact.DISPLAY_FIELDS = [].concat(
	ZmContact.ADDRESS_FIELDS,
	ZmContact.EMAIL_FIELDS,
	ZmContact.IM_FIELDS,
	ZmContact.OTHER_FIELDS,
	ZmContact.PHONE_FIELDS,
	ZmContact.PRIMARY_FIELDS,
	ZmContact.URL_FIELDS
);

ZmContact.IGNORE_FIELDS = [].concat(
	ZmContact.GAL_FIELDS,
	ZmContact.MYCARD_FIELDS,
	ZmContact.X_FIELDS,
	[ZmContact.F_imagepart]
);

ZmContact.ALL_FIELDS = [].concat(
	ZmContact.DISPLAY_FIELDS, ZmContact.IGNORE_FIELDS
);

ZmContact.IS_DATE = {};
ZmContact.IS_DATE[ZmContact.F_birthday] = true;
ZmContact.IS_DATE[ZmContact.F_anniversary] = true;

ZmContact.IS_IGNORE = AjxUtil.arrayAsHash(ZmContact.IGNORE_FIELDS);

// number of distribution list members to fetch at a time
ZmContact.DL_PAGE_SIZE = 100;

ZmContact.GROUP_CONTACT_REF = "C";
ZmContact.GROUP_GAL_REF = "G";
ZmContact.GROUP_INLINE_REF = "I";	
}; // updateFieldConstants()
ZmContact.updateFieldConstants();

/**
 * This structure can be queried to determine if the first
 * entry in a multi-value entry is suffixed with "1". Most
 * attributes add a numerical suffix to all but the first
 * entry.
 * <p>
 * <strong>Note:</strong>
 * In most cases, {@link ZmContact#getAttributeName} is a better choice.
 */
ZmContact.IS_ADDONE = {};
ZmContact.IS_ADDONE[ZmContact.F_custom] = true;
ZmContact.IS_ADDONE[ZmContact.F_imAddress] = true;
ZmContact.IS_ADDONE[ZmContact.X_outlookUserField] = true;

/**
 * Gets an indexed attribute name taking into account if the field
 * with index 1 should append the "1" or not. Code should call this
 * function in lieu of accessing {@link ZmContact.IS_ADDONE} directly.
 */
ZmContact.getAttributeName = function(name, index) {
	index = index || 1;
	return index > 1 || ZmContact.IS_ADDONE[name] ? name+index : name;
};

/**
 * Returns a string representation of the object.
 * 
 * @return		{String}		a string representation of the object
 */
ZmContact.prototype.toString =
function() {
	return "ZmContact";
};

// Class methods

/**
 * Creates a contact from an XML node.
 *
 * @param {Object}	node		a "cn" XML node
 * @param {Hash}	args		args to pass to the constructor
 * @return	{ZmContact}	the contact
 */
ZmContact.createFromDom =
function(node, args) {
	// check global cache for this item first
	var contact = appCtxt.cacheGet(node.id);

	// make sure the revision hasnt changed, otherwise contact is out of date
	if (contact == null || (contact && contact.rev != node.rev)) {
		contact = new ZmContact(node.id, args.list);
		if (args.isGal) {
			contact.isGal = args.isGal;
		}
		contact._loadFromDom(node);
		//update the canonical list
		appCtxt.getApp(ZmApp.CONTACTS).getContactList().add(contact);
	} else {
		if (node.m) {
			contact.attr[ZmContact.F_groups] = node.m;
		}
		if (node.ref) {
			contact.ref = node.ref;
		}
		if (node.tn) {
			contact._parseTagNames(node.tn);
		}
		AjxUtil.hashUpdate(contact.attr, node._attrs);	// merge new attrs just in case we don't have them
		contact.list = args.list || new ZmContactList(null);
		contact._list = {};
		contact._list[contact.list.id] = true;
	}
	//S/MIME: If user certificate is present, include it into contact's object
	if (node.certificate) {
		contact.certificate = node.certificate;
	}

	return contact;
};

/**
 * Compares two contacts based on how they are filed. Intended for use by
 * sort methods.
 *
 * @param {ZmContact}		a		a contact
 * @param {ZmContact}		b		a contact
 * @return	{int}	0 if the contacts are the same; 1 if "a" is before "b"; -1 if "b" is before "a"
 */
ZmContact.compareByFileAs =
function(a, b) {
	var aFileAs = (a instanceof ZmContact) ? a.getFileAs(true) : ZmContact.computeFileAs(a._attrs).toLowerCase();
	var bFileAs = (b instanceof ZmContact) ? b.getFileAs(true) : ZmContact.computeFileAs(b._attrs).toLowerCase();

	if (!bFileAs || (aFileAs > bFileAs)) return 1;
	if (aFileAs < bFileAs) return -1;
	return 0;
};

/**
 * Figures out the filing string for the contact according to the chosen method.
 *
 * @param {ZmContact|Hash}	contact		a contact or a hash of contact attributes
 */
ZmContact.computeFileAs =
function(contact) {
	/*
	 * Bug 98176: To keep the same logic of generating the FileAs contact
	 *    label string between the Ajax client, and HTML client, when the
	 *    computeFileAs(), and fileAs*() functions are modified, please
	 *    change the corresponding functions defined in the autoComplete.tag
	 */
	var attr = (contact instanceof ZmContact) ? contact.getAttrs() : contact;
	if (!attr) return;

	if (attr[ZmContact.F_dlDisplayName]) {
		//this is only DL case. But since this is sometimes just the attrs,
		//I can't always use isDistributionList method.
		return attr[ZmContact.F_dlDisplayName];
	}

	var val = parseInt(attr.fileAs);
	var fa;
	var idx = 0;

	switch (val) {
		case ZmContact.FA_LAST_C_FIRST: 										// Last, First
		default: {
			// if GAL contact, use full name instead (bug fix #4850,4009)
			if (contact && contact.isGal) {
				if (attr.fullName) { // bug fix #27428 - if fullName is Array, return first
					return (attr.fullName instanceof Array) ? attr.fullName[0] : attr.fullName;
				}
				return ((attr.email instanceof Array) ? attr.email[0] : attr.email);
			}
			fa = ZmContact.fileAsLastFirst(attr.firstName, attr.lastName, attr.fullName, attr.nickname);
		}
		break;

		case ZmContact.FA_FIRST_LAST: { 										// First Last
			fa = ZmContact.fileAsFirstLast(attr.firstName, attr.lastName, attr.fullName, attr.nickname);
		}
		break;

		case ZmContact.FA_COMPANY: {											// Company
			if (attr.company) fa = attr.company;
		}
		break;

		case ZmContact.FA_LAST_C_FIRST_COMPANY: {								// Last, First (Company)
			var name = ZmContact.fileAsLastFirst(attr.firstName, attr.lastName, attr.fullName, attr.nickname);
			fa = ZmContact.fileAsNameCompany(name, attr.company);
		}
		break;

		case ZmContact.FA_FIRST_LAST_COMPANY: {									// First Last (Company)
			var name = ZmContact.fileAsFirstLast(attr.firstName, attr.lastName, attr.fullName, attr.nickname);
			fa = ZmContact.fileAsNameCompany(name, attr.company);
		}
		break;

		case ZmContact.FA_COMPANY_LAST_C_FIRST: {								// Company (Last, First)
			var name = ZmContact.fileAsLastFirst(attr.firstName, attr.lastName);
			fa = ZmContact.fileAsCompanyName(name, attr.company);
		}
		break;

		case ZmContact.FA_COMPANY_FIRST_LAST: {									// Company (First Last)
			var name = ZmContact.fileAsFirstLast(attr.firstName, attr.lastName);
			fa = ZmContact.fileAsCompanyName(name, attr.company);
		}
		break;

		case ZmContact.FA_CUSTOM: {												// custom looks like this: "8:foobar"
			return attr.fileAs.substring(2);
		}
		break;
	}
	return fa || attr.fullName || "";
};

/**
 * Name printing helper "First Last".
 * 
 * @param	{String}	first		the first name
 * @param	{String}	last		the last name
 * @param	{String}	fullname		the fullname
 * @param	{String}	nickname		the nickname
 * @return	{String}	the name format
 */
ZmContact.fileAsFirstLast =
function(first, last, fullname, nickname) {
	if (first && last)
		return AjxMessageFormat.format(ZmMsg.fileAsFirstLast, [first, last]);
	return first || last || fullname || nickname || "";
};

/**
 * Name printing helper "Last, First".
 * 
 * @param	{String}	first		the first name
 * @param	{String}	last		the last name
 * @param	{String}	fullname		the fullname
 * @param	{String}	nickname		the nickname
 * @return	{String}	the name format
 */
ZmContact.fileAsLastFirst =
function(first, last, fullname, nickname) {
	if (first && last)
		return AjxMessageFormat.format(ZmMsg.fileAsLastFirst, [first, last]);
	return last || first || fullname || nickname || "";
};

/**
 * Name printing helper "Name (Company)".
 *
 * @param	{String}	name		the contact name
 * @param	{String}	company		the company
 * @return	{String}	the name format
 */
ZmContact.fileAsNameCompany =
function(name, company) {
	if (name && company)
		return AjxMessageFormat.format(ZmMsg.fileAsNameCompany, [name, company]);
	if (company)
		return AjxMessageFormat.format(ZmMsg.fileAsCompanyAsSecondaryOnly, [company]);
	return name;
};

/**
 * Name printing helper "Company (Name)".
 * 
 * @param	{String}	name		the contact name
 * @param	{String}	company		the company
 * @return	{String}	the name format
 */
ZmContact.fileAsCompanyName =
function(name, company) {
	if (company && name)
		return AjxMessageFormat.format(ZmMsg.fileAsCompanyName, [name, company]);
	if (name)
		return AjxMessageFormat.format(ZmMsg.fileAsNameAsSecondaryOnly, [name]);
	return company;
};

/**
 * Computes the custom file as string by prepending "8:" to the given custom fileAs string.
 * 
 * @param {Hash}	customFileAs	a set of contact attributes
 * @return	{String}	the name format
 */
ZmContact.computeCustomFileAs =
function(customFileAs) {
	return [ZmContact.FA_CUSTOM, ":", customFileAs].join("");
};

/*
 * 
 * These next few static methods handle a contact that is either an anonymous
 * object or an actual ZmContact. The former is used to optimize loading. The
 * anonymous object is upgraded to a ZmContact when needed.
 *  
 */

/**
 * Gets an attribute.
 * 
 * @param	{ZmContact}	contact		the contact
 * @param	{String}	attr		the attribute
 * @return	{Object}	the attribute value or <code>null</code> for none
 */
ZmContact.getAttr =
function(contact, attr) {
	return (contact instanceof ZmContact)
		? contact.getAttr(attr)
		: (contact && contact._attrs) ? contact._attrs[attr] : null;
};

/**
 * returns the prefix of a string in the format "abc123". (would return "abc"). If the string is all number, it's a special case and returns the string itself. e.g. "234" would return "234".
 */
ZmContact.getPrefix = function(s) {
	var trimmed = s.replace(/\d+$/, "");
	if (trimmed === "") {
		//number only - don't trim. The number is the prefix.
		return s;
	}
	return trimmed;
};

/**
 * Normalizes the numbering of the given attribute names and
 * returns a new object with the re-numbered attributes. For
 * example, if the attributes contains a "foo2" but no "foo",
 * then the "foo2" attribute will be renamed to "foo" in the
 * returned object.
 *
 * @param {Hash}	attrs  a hash of attributes to normalize.
 * @param {String}	[prefix] if specified, only the the attributes that match the given prefix will be returned
 * @param {Array}	[ignore] if specified, the attributes that are present in the array will not be normalized
 * @return	{Hash}	a hash of normalized attributes
 */
ZmContact.getNormalizedAttrs = function(attrs, prefix, ignore) {
	var nattrs = {};
	if (attrs) {
		// normalize attribute numbering
		var names = AjxUtil.keys(attrs);
		names.sort(ZmContact.__BY_ATTRIBUTE);
		var a = {};
		for (var i = 0; i < names.length; i++) {
			var name = names[i];
			// get current count
			var nprefix = ZmContact.getPrefix(name);
			if (prefix && prefix != nprefix) continue;
			if (AjxUtil.isArray(ignore) && AjxUtil.indexOf(ignore, nprefix)!=-1) {
				nattrs[name] = attrs[name];
			} else {
				if (!a[nprefix]) a[nprefix] = 0;
				// normalize, if needed
				var nname = ZmContact.getAttributeName(nprefix, ++a[nprefix]);
				nattrs[nname] = attrs[name];
			}
		}
	}
	return nattrs;
};

ZmContact.__RE_ATTRIBUTE = /^(.*?)(\d+)$/;
ZmContact.__BY_ATTRIBUTE = function(a, b) {
	var aa = a.match(ZmContact.__RE_ATTRIBUTE) || [a,a,1];
	var bb = b.match(ZmContact.__RE_ATTRIBUTE) || [b,b,1];
	return aa[1] == bb[1] ? Number(aa[2]) - Number(bb[2]) : aa[1].localeCompare(bb[1]);
};

/**
 * Sets the attribute.
 * 
 * @param	{ZmContact}	contact		the contact
 * @param	{String}	attr		the attribute
 * @param	{Object}	value		the attribute value
 */
ZmContact.setAttr =
function(contact, attr, value) {
	if (contact instanceof ZmContact)
		contact.setAttr(attr, value);
	else
		contact._attrs[attr] = value;
};

/**
 * Checks if the contact is in the trash.
 * 
 * @param	{ZmContact}	contact		the contact
 * @return	{Boolean}	<code>true</code> if in trash
 */
ZmContact.isInTrash =
function(contact) {
	var folderId = (contact instanceof ZmContact) ? contact.folderId : contact.l;
	var folder = appCtxt.getById(folderId);
	return (folder && folder.isInTrash());
};

/**
 * @private
 */
ZmContact.prototype.load =
function(callback, errorCallback, batchCmd, deref) {
	var jsonObj = {GetContactsRequest:{_jsns:"urn:zimbraMail"}};
	if (deref) {
		jsonObj.GetContactsRequest.derefGroupMember = "1";
	}
	jsonObj.GetContactsRequest.returnCertInfo = "1"; //Fix for: ZCS-999 and ZCS-991
	var request = jsonObj.GetContactsRequest;
	request.cn = [{id:this.id}];

	var respCallback = new AjxCallback(this, this._handleLoadResponse, [callback]);

	if (batchCmd) {
		var jsonObj = {GetContactsRequest:{_jsns:"urn:zimbraMail"}};
		if (deref) {
			jsonObj.GetContactsRequest.derefGroupMember = "1";
		}
		jsonObj.GetContactsRequest.returnCertInfo = "1"; //Fix for: ZCS-999 and ZCS-991
		jsonObj.GetContactsRequest.cn = {id:this.id};
		batchCmd.addRequestParams(jsonObj, respCallback, errorCallback);
	} else {
		appCtxt.getAppController().sendRequest({jsonObj:jsonObj,
												asyncMode:true,
												callback:respCallback,
												errorCallback:errorCallback});
	}
};

/**
 * @private
 */
ZmContact.prototype._handleLoadResponse =
function(callback, result) {
	var resp = result.getResponse().GetContactsResponse;

	// for now, we just assume only one contact was requested at a time
	var contact = resp.cn[0];
	this.attr = contact._attrs;
	if (contact.m) {
		for (var i = 0; i < contact.m.length; i++) {
			//cache contacts from contact groups (e.g. GAL contacts, shared contacts have not already been cached)
			var member = contact.m[i];
			var isGal = false;
			if (member.type == ZmContact.GROUP_GAL_REF) {
				isGal = true;
			}
			if (member.cn && member.cn.length > 0) {
				var memberContact = member.cn[0];
				memberContact.ref = memberContact.ref || (isGal && member.value); //we sometimes don't get "ref" but the "value" for GAL is the ref.
				var loadMember = ZmContact.createFromDom(memberContact, {list: this.list, isGal: isGal}); //pass GAL so fileAS gets set correctly
				loadMember.isDL = isGal && loadMember.attr[ZmContact.F_type] == "group";
				appCtxt.cacheSet(member.value, loadMember);
			}
			
		}
		this._loadFromDom(contact); //load group
	}
	this.isLoaded = true;
	if (callback) {
		callback.run(contact, this);
	}
};

/**
 * @private
 */
ZmContact.prototype.clear =
function() {
	// bug fix #41666 - override base class method and do nothing
};

/**
 * Checks if the contact attributes are empty.
 * 
 * @return	{Boolean}	<code>true</code> if empty
 */
ZmContact.prototype.isEmpty =
function() {
	for (var i in this.attr) {
		return false;
	}
	return true;
};

/**
 * Checks if the contact is shared.
 * 
 * @return	{Boolean}	<code>true</code> if shared
 */
ZmContact.prototype.isShared =
function() {
	return this.addrbook && this.addrbook.link;
};

/**
 * Checks if the contact is read-only.
 * 
 * @return	{Boolean}	<code>true</code> if read-only
 */
ZmContact.prototype.isReadOnly =
function() {
	if (this.isGal) { return true; }

	return this.isShared()
		? this.addrbook && this.addrbook.isReadOnly()
		: false;
};

/**
 * Checks if the contact is locked. This is different for DLs than read-only.
 *
 * @return	{Boolean}	<code>true</code> if read-only
 */
ZmContact.prototype.isLocked =
function() {
	if (!this.isDistributionList()) {
		return this.isReadOnly();
	}
	if (!this.dlInfo) {
		return false; //rare case after editing by an owner if the fileAsChanged, the new dl Info still not read, and the layout re-done. So don't show the lock.
	}
	var dlInfo = this.dlInfo;
	if (dlInfo.isOwner) {
		return false;
	}
	if (dlInfo.isMember) {
    	return dlInfo.unsubscriptionPolicy == ZmContactSplitView.SUBSCRIPTION_POLICY_REJECT;
	}
	return dlInfo.subscriptionPolicy == ZmContactSplitView.SUBSCRIPTION_POLICY_REJECT;
};

/**
 * Checks if the contact is a group.
 * 
 * @return	{Boolean}	<code>true</code> if a group
 */
ZmContact.prototype.isGroup =
function() {
	return this.getAttr(ZmContact.F_type) == "group" || this.type == ZmItem.GROUP;
};

/**
 * Checks if the contact is a DL.
 *
 * @return	{Boolean}	<code>true</code> if a group
 */
ZmContact.prototype.isDistributionList =
function() {
	return this.isGal && this.isGroup();
};


// parses "groups" attr into AjxEmailAddress objects stored in 3 vectors (all, good, and bad)
/**
 * Gets the group members.
 *
 * @return	{AjxVector}		the group members or <code>null</code> if not group
 */
ZmContact.prototype.getGroupMembers =
function() {
	var allMembers = this.getAllGroupMembers();
	var addrs = [];
	for (var i = 0; i < allMembers.length; i++) {
		addrs.push(allMembers[i].toString());
	}
	return AjxEmailAddress.parseEmailString(addrs.join(", "));
};	

/**
 * parses "groups" attr into an AjxEmailAddress with a few extra attributes (see ZmContactsHelper._wrapInlineContact)
 * 
 * @return	{AjxVector}		the group members or <code>null</code> if not group
 */
ZmContact.prototype.getAllGroupMembers =
function() {

	if (this.isDistributionList()) {
		return this.dlMembers;
	}

	var addrs = [];

	var groupMembers = this.attr[ZmContact.F_groups];
	if (!groupMembers){
		return AjxEmailAddress.parseEmailString(this.attr[ZmContact.F_email]);  //I doubt this is needed or works correctly, but I keep this logic from before. If we don't have the group members, how can we return the group email instead?
	}
	for (var i = 0; i < groupMembers.length; i++) {
		var member = groupMembers[i];
		var type = member.type;
		var value = member.value;
		if (type == ZmContact.GROUP_INLINE_REF) {
			addrs.push(ZmContactsHelper._wrapInlineContact(value));
		}
		else {
			var contact = ZmContact.getContactFromCache(value);	 //TODO: handle contacts not cached?
			if (!contact) {
				DBG.println(AjxDebug.DBG1, "Disregarding uncached contact: " + value);
				continue;
			}
			var ajxEmailAddress = ZmContactsHelper._wrapContact(contact);
			if (ajxEmailAddress && type === ZmContact.GROUP_CONTACT_REF) {
				ajxEmailAddress.groupRefValue = value; //don't normalize value
			}
			if (ajxEmailAddress) {
				addrs.push(ajxEmailAddress);
			}
		}
	}
	return addrs;
};


ZmContact.prototype.gatherExtraDlStuff =
function(callback) {
	if (this.dlInfo && !this.dlInfo.isMinimal) {
		//already there, skip to next step, loading DL Members
		this.loadDlMembers(callback);
		return;
	}
	var callbackFromGettingInfo = this._handleGetDlInfoResponse.bind(this, callback);
	this.loadDlInfo(callbackFromGettingInfo);
};


ZmContact.prototype._handleGetDlInfoResponse =
function(callback, result) {
	var response = result._data.GetDistributionListResponse;
	var dl = response.dl[0];
	var attrs = dl._attrs;
	var isMember = dl.isMember;
	var isOwner = dl.isOwner;
	var mailPolicySpecificMailers = [];
	this.dlInfo = {	isMember: isMember,
						isOwner: isOwner,
						subscriptionPolicy: attrs.zimbraDistributionListSubscriptionPolicy,
						unsubscriptionPolicy: attrs.zimbraDistributionListUnsubscriptionPolicy,
						description: attrs.description || "",
						displayName: attrs.displayName || "",
						notes: attrs.zimbraNotes || "",
						hideInGal: attrs.zimbraHideInGal == "TRUE",
						mailPolicy: isOwner && this._getMailPolicy(dl, mailPolicySpecificMailers),
						owners: isOwner && this._getOwners(dl)};
	this.dlInfo.mailPolicySpecificMailers = mailPolicySpecificMailers;

	this.loadDlMembers(callback);
};

ZmContact.prototype.loadDlMembers =
function(callback) {
	if ((!appCtxt.get("EXPAND_DL_ENABLED") || this.dlInfo.hideInGal) && !this.dlInfo.isOwner) {
		// can't get members if dl has zimbraHideInGal true, and not owner
		//also, if zimbraFeatureDistributionListExpandMembersEnabled is false - also do not show the members (again unless it's the owner)
		this.dlMembers = [];
		if (callback) {
			callback();
		}
		return;
	}
	if (this.dlMembers) {
		//already there - just callback
		if (callback) {
			callback();
		}
		return;
	}
	var respCallback = this._handleGetDlMembersResponse.bind(this, callback);
	this.getAllDLMembers(respCallback);
};


ZmContact.prototype._handleGetDlMembersResponse =
function(callback, result) {
	var list = result.list;
	if (!list) {
		this.dlMembers = [];
		callback();
		return;
	}
	var members = [];
	for (var i = 0; i < list.length; i++) {
		members.push({type: ZmContact.GROUP_INLINE_REF,
						value: list[i],
						address: list[i]});
	}

	this.dlMembers = members;
	callback();
};

ZmContact.prototype._getOwners =
function(dl) {
	var owners = dl.owners[0].owner;
	var ownersArray = [];
	for (var i = 0; i < owners.length; i++) {
		var owner = owners[i].name;
		ownersArray.push(owner); //just the email address, I think and hope.
	}
	return ownersArray;
};

ZmContact.prototype._getMailPolicy =
function(dl, specificMailers) {
	var mailPolicy;

	var rights = dl.rights[0].right;
	var right = rights[0];
	var grantees = right.grantee;
	if (!grantees) {
		return ZmGroupView.MAIL_POLICY_ANYONE;
	}
	for (var i = 0; i < grantees.length; i++) {
		var grantee = grantees[i];

		mailPolicy = ZmGroupView.GRANTEE_TYPE_TO_MAIL_POLICY_MAP[grantee.type];

		if (mailPolicy == ZmGroupView.MAIL_POLICY_SPECIFIC) {
			specificMailers.push(grantee.name);
		}
		else if (mailPolicy == ZmGroupView.MAIL_POLICY_ANYONE) {
			break;
		}
		else if (mailPolicy == ZmGroupView.MAIL_POLICY_INTERNAL) {
			break;
		}
		else if (mailPolicy == ZmGroupView.MAIL_POLICY_MEMBERS) {
			if (grantee.name == this.getEmail()) {
				//this means only members of this DL can send.
				break;
			}
			else {
				//must be another DL, and we do allow it, so treat it as regular user.
				specificMailers.push(grantee.name);
				mailPolicy = ZmGroupView.MAIL_POLICY_SPECIFIC;
			}
		}
	}
	mailPolicy = mailPolicy || ZmGroupView.MAIL_POLICY_ANYONE;

	return mailPolicy;
};


ZmContact.prototype.loadDlInfo =
function(callback) {
	var soapDoc = AjxSoapDoc.create("GetDistributionListRequest", "urn:zimbraAccount", null);
	soapDoc.setMethodAttribute("needOwners", "1");
	soapDoc.setMethodAttribute("needRights", "sendToDistList");
	var elBy = soapDoc.set("dl", this.getEmail());
	elBy.setAttribute("by", "name");

	appCtxt.getAppController().sendRequest({soapDoc: soapDoc, asyncMode: true, callback: callback});
};

ZmContact.prototype.toggleSubscription =
function(callback) {
	var soapDoc = AjxSoapDoc.create("SubscribeDistributionListRequest", "urn:zimbraAccount", null);
	soapDoc.setMethodAttribute("op", this.dlInfo.isMember ? "unsubscribe" : "subscribe");
	var elBy = soapDoc.set("dl", this.getEmail());
	elBy.setAttribute("by", "name");
	appCtxt.getAppController().sendRequest({soapDoc: soapDoc, asyncMode: true, callback: callback});
};



/**
 *  Returns the contact id.  If includeUserZid is true it will return the format zid:id
 * @param includeUserZid {boolean} true to include the zid prefix for the contact id
 * @return {String} contact id string
 */
ZmContact.prototype.getId = 
function(includeUserZid) {

	if (includeUserZid) {
		return this.isShared() ? this.id : appCtxt.accountList.mainAccount.id + ":" + this.id; 
	}
	
	return this.id;
};
/**
 * Gets the icon.
 * @param 	{ZmAddrBook} addrBook	address book of contact 
 * @return	{String}	the icon
 */
ZmContact.prototype.getIcon =
function(addrBook) {
	if (this.isDistributionList()) 						{ return "DistributionList"; }
	if (this.isGal)										{ return "GALContact"; }
	if (this.isShared() || (addrBook && addrBook.link))	{ return "SharedContact"; }
	if (this.isGroup())									{ return "Group"; }
	return "Contact";
};

ZmContact.prototype.getIconLarge =
function() {
	if (this.isDistributionList()) {
		return "Group_48";
	}
	//todo - get a big version of ImgGalContact.png
//	if (this.isGal) {
//	}
	return "Person_48";
};

/**
 * Gets the folder id.
 * 
 * @return	{String}		the folder id	
 */
ZmContact.prototype.getFolderId =
function() {
	return this.isShared()
		? this.folderId.split(":")[0]
		: this.folderId;
};

/**
 * Gets the attribute.
 * 
 * @param	{String}	name		the attribute name
 * @return	{String}	the value
 */
ZmContact.prototype.getAttr =
function(name) {
	var val = this.attr[name];
	return val ? ((val instanceof Array) ? val[0] : val) : "";
};

/**
 * Sets the attribute.
 * 
 * @param	{String}	name		the attribute name
 * @param	{String}	value		the attribute value
 */
ZmContact.prototype.setAttr =
function(name, value) {
	this.attr[name] = value;
};

/**
 * Sets the participant status.
 *
 * @param	{String}	value the participant status value
 */
ZmContact.prototype.setParticipantStatus =
function(ptst) {
	this.participantStatus = ptst;
};

/**
 * gets the participant status.
 *
 * @return	{String}    the value
 */
ZmContact.prototype.getParticipantStatus =
function() {
	return this.participantStatus;
};

/**
 * Sets the participant role.
 *
 * @param	{String}	value the participant role value
 */
ZmContact.prototype.setParticipantRole =
function(role) {
	this.participantRole = role;
};

/**
 * gets the participant role.
 *
 * @return	{String}    the value
 */
ZmContact.prototype.getParticipantRole =
function() {
	return this.participantRole;
};

/**
 * Removes the attribute.
 * 
 * @param	{String}	name		the attribute name
 */
ZmContact.prototype.removeAttr =
function(name) {
	delete this.attr[name];
};

/**
 * Gets the contact attributes.
 *
 * @param {String}	[prefix] if specified, only the the attributes that match the given prefix will be returned
 * @return	{Hash}	a hash of attribute/value pairs
 */
ZmContact.prototype.getAttrs = function(prefix) {
	var attrs = this.attr;
	if (prefix) {
		attrs = {};
		for (var aname in this.attr) {
			var namePrefix = ZmContact.getPrefix(aname);
			if (namePrefix === prefix) {
				attrs[aname] = this.attr[aname];
			}
		}
	}
	return attrs;
};

/**
 * Gets a normalized set of attributes where the attribute
 * names have been re-numbered as needed. For example, if the
 * attributes contains a "foo2" but no "foo", then the "foo2"
 * attribute will be renamed to "foo" in the returned object.
 * <p>
 * <strong>Note:</strong>
 * This method is expensive so should be called once and
 * cached temporarily as needed instead of being called
 * for each normalized attribute that is needed.
 * 
 * @param {String}	[prefix]		if specified, only the
 *                        the attributes that match the given
 *                        prefix will be returned.
 * @return	{Hash}	a hash of attribute/value pairs
 */
ZmContact.prototype.getNormalizedAttrs = function(prefix) {
	return ZmContact.getNormalizedAttrs(this.attr, prefix, ZmContact.IGNORE_NORMALIZATION);
};

/**
* Creates a contact from the given set of attributes. Used to create contacts on
* the fly (rather than by loading them). This method is called by a list's <code>create()</code>
* method.
* <p>
* If this is a GAL contact, we assume it is being added to the contact list.</p>
*
* @param {Hash}	attr			the attribute/value pairs for this contact
* @param {ZmBatchCommand}	batchCmd	the batch command that contains this request
* @param {boolean} isAutoCreate true if this is a auto create and toast message should not be shown
*/
ZmContact.prototype.create =
function(attr, batchCmd, isAutoCreate) {

	if (this.isDistributionList()) {
		this._createDl(attr);
		return;
	}

	var jsonObj = {CreateContactRequest:{_jsns:"urn:zimbraMail"}};
	var request = jsonObj.CreateContactRequest;
	var cn = request.cn = {};

	var folderId = attr[ZmContact.F_folderId] || ZmFolder.ID_CONTACTS;
	var folder = appCtxt.getById(folderId);
	if (folder && folder.isRemote()) {
		folderId = folder.getRemoteId();
	}
	cn.l = folderId;
	cn.a = [];
	cn.m = [];

	for (var name in attr) {
		if (name == ZmContact.F_folderId ||
			name == "objectClass" ||
			name == "zimbraId" ||
			name == "createTimeStamp" ||
			name == "modifyTimeStamp") { continue; }

		if (name == ZmContact.F_groups) {
			this._addContactGroupAttr(cn, attr);
		}
		else {
			this._addRequestAttr(cn, name, attr[name]);
		}
	}

	this._addRequestAttr(cn, ZmContact.X_fullName, ZmContact.computeFileAs(attr));

	var respCallback = new AjxCallback(this, this._handleResponseCreate, [attr, batchCmd != null, isAutoCreate]);

	if (batchCmd) {
		batchCmd.addRequestParams(jsonObj, respCallback);
	} else {
		appCtxt.getAppController().sendRequest({jsonObj:jsonObj, asyncMode:true, callback:respCallback});
	}
};

/**
 * @private
 */
ZmContact.prototype._handleResponseCreate =
function(attr, isBatchMode, isAutoCreate, result) {
	// dont bother processing creates when in batch mode (just let create
	// notifications handle them)
	if (isBatchMode) { return; }

	var resp = result.getResponse().CreateContactResponse;
	cn = resp ? resp.cn[0] : null;
	var id = cn ? cn.id : null;
	if (id) {
		this._fileAs = null;
		this._fullName = null;
		this.id = id;
		this.modified = cn.md;
		this.folderId = cn.l || ZmOrganizer.ID_ADDRBOOK;
		for (var a in attr) {
			if (!(attr[a] == undefined || attr[a] == ''))
				this.setAttr(a, attr[a]);
		}
		var groupMembers = cn ? cn.m : null;
		if (groupMembers) {
			this.attr[ZmContact.F_groups] = groupMembers;
			cn._attrs[ZmContact.F_groups] = groupMembers;
		}
		if (!isAutoCreate) {
			var msg = this.isGroup() ? ZmMsg.groupCreated : ZmMsg.contactCreated;
			appCtxt.getAppController().setStatusMsg(msg);
		}
		//update the canonical list. (this includes adding to the _idHash like before (bug 44132) calling updateIdHash. But calling that left the list inconcistant.
		appCtxt.getApp(ZmApp.CONTACTS).getContactList().add(cn);
	} else {
		var msg = this.isGroup() ? ZmMsg.errorCreateGroup : ZmMsg.errorCreateContact;
		var detail = ZmMsg.errorTryAgain + "\n" + ZmMsg.errorContact;
		appCtxt.getAppController().setStatusMsg(msg, ZmStatusView.LEVEL_CRITICAL, detail);
	}
};

/**
 * Creates a contct from a VCF part of a message.
 * 
 * @param	{String}	msgId		the message
 * @param	{String}	vcardPartId	the vcard part id
 */
ZmContact.prototype.createFromVCard =
function(msgId, vcardPartId) {
	var jsonObj = {CreateContactRequest:{_jsns:"urn:zimbraMail"}};
	var cn = jsonObj.CreateContactRequest.cn = {l:ZmFolder.ID_CONTACTS};
	cn.vcard = {mid:msgId, part:vcardPartId};

	var params = {
		jsonObj: jsonObj,
		asyncMode: true,
		callback: (new AjxCallback(this, this._handleResponseCreateVCard)),
		errorCallback: (new AjxCallback(this, this._handleErrorCreateVCard))
	};

	appCtxt.getAppController().sendRequest(params);
};

/**
 * @private
 */
ZmContact.prototype._handleResponseCreateVCard =
function(result) {
	appCtxt.getAppController().setStatusMsg(ZmMsg.contactCreated);
};

/**
 * @private
 */
ZmContact.prototype._handleErrorCreateVCard =
function(ex) {
	appCtxt.getAppController().setStatusMsg(ZmMsg.errorCreateContact, ZmStatusView.LEVEL_CRITICAL);
};

/**
 * Updates contact attributes.
 *
 * @param {Hash}	attr		a set of attributes and new values
 * @param {AjxCallback}	callback	the callback
 * @param {boolean} isAutoSave  true if it is a auto save and toast should not be displayed.
 */
ZmContact.prototype.modify =
function(attr, callback, isAutoSave, batchCmd) {
	if (this.isDistributionList()) {
		this._modifyDl(attr);
		return;
	}
	if (this.list.isGal) { return; }

	// change force to 0 and put up dialog if we get a MODIFY_CONFLICT fault?
	var jsonObj = {ModifyContactRequest:{_jsns:"urn:zimbraMail", replace:"0", force:"1"}};
	var cn = jsonObj.ModifyContactRequest.cn = {id:this.id};
	cn.a = [];
	cn.m = [];
	var continueRequest = false;
	
	for (var name in attr) {
		if (name == ZmContact.F_folderId) { continue; }
		if (name == ZmContact.F_groups) {
			this._addContactGroupAttr(cn, attr);	
		}
		else {
			this._addRequestAttr(cn, name, (attr[name] && attr[name].value) || attr[name]);
		}
		continueRequest = true;
	}

    // bug: 45026
    if (ZmContact.F_firstName in attr || ZmContact.F_lastName in attr || ZmContact.F_company in attr || ZmContact.X_fileAs in attr) {
        var contact = {};
        var fields = [ZmContact.F_firstName, ZmContact.F_lastName, ZmContact.F_company, ZmContact.X_fileAs];
        for (var i = 0; i < fields.length; i++) {
            var field = fields[i];
            var value = attr[field];
            contact[field] = value != null ? value : this.getAttr(field);
        }
        var fullName = ZmContact.computeFileAs(contact); 
        this._addRequestAttr(cn, ZmContact.X_fullName, fullName);
    }

	if (continueRequest) {
		if (batchCmd) {
			batchCmd.addRequestParams(jsonObj, null, null); //no need for response callback for current use-case (batch modifying zimlet image)
		}
		else {
			var respCallback = this._handleResponseModify.bind(this, attr, callback, isAutoSave);
			appCtxt.getAppController().sendRequest({jsonObj: jsonObj, asyncMode: true, callback: respCallback});
		}

	} else {
		if (attr[ZmContact.F_folderId]) {
			this._setFolder(attr[ZmContact.F_folderId]);
		}
	}
};

ZmContact.prototype._createDl =
function(attr) {

	this.attr = attr; //this is mainly important for the email. attr is not set before this.

	var createDlReq = this._getCreateDlReq(attr);

	var reqs = [];

	this._addMemberModsReqs(reqs, attr);

	this._addMailPolicyAndOwnersReqs(reqs, attr);

	var jsonObj = {
		BatchRequest: {
			_jsns: "urn:zimbra",
			CreateDistributionListRequest: createDlReq,
			DistributionListActionRequest: reqs
		}
	};
	var respCallback = this._createDlResponseHandler.bind(this);
	appCtxt.getAppController().sendRequest({jsonObj: jsonObj, asyncMode: true, callback: respCallback});
	
};

ZmContact.prototype._addMailPolicyAndOwnersReqs =
function(reqs, attr) {

	var mailPolicy = attr[ZmContact.F_dlMailPolicy];
	if (mailPolicy) {
		reqs.push(this._getSetMailPolicyReq(mailPolicy, attr[ZmContact.F_dlMailPolicySpecificMailers]));
	}

	var listOwners = attr[ZmContact.F_dlListOwners];
	if (listOwners) {
		reqs.push(this._getSetOwnersReq(listOwners));
	}


};



ZmContact.prototype._addMemberModsReqs =
function(reqs, attr) {
	var memberModifications = attr[ZmContact.F_groups];
	var adds = [];
	var removes = [];
	if (memberModifications) {
		for (var i = 0; i < memberModifications.length; i++) {
			var mod = memberModifications[i];
			var col = (mod.op == "+" ? adds : removes);
			col.push(mod);
		}
	}

	if (adds.length > 0) {
		reqs.push(this._getAddOrRemoveReq(adds, true));
	}
	if (removes.length > 0) {
		reqs.push(this._getAddOrRemoveReq(removes, false));
	}
};

ZmContact.prototype._modifyDl =
function(attr) {
	var reqs = [];

	var newEmail = attr[ZmContact.F_email];

	var emailChanged = false;
	if (newEmail !== undefined) {
		emailChanged = true;
		reqs.push(this._getRenameDlReq(newEmail));
		this.setAttr(ZmContact.F_email, newEmail);
	}

	var modDlReq = this._getModifyDlAttributesReq(attr);
	if (modDlReq) {
		reqs.push(modDlReq);
	}

	var displayName = attr[ZmContact.F_dlDisplayName];
	if (displayName !== undefined) {
		this.setAttr(ZmContact.F_dlDisplayName, displayName);
	}

	var oldFileAs = this.getFileAs();
	this._resetCachedFields();
	var fileAsChanged = oldFileAs != this.getFileAs();

	this._addMemberModsReqs(reqs, attr);

	this._addMailPolicyAndOwnersReqs(reqs, attr);

	if (reqs.length == 0) {
		this._modifyDlResponseHandler(false, null); //pretend it was saved
		return;
	}
	var jsonObj = {
		BatchRequest: {
			_jsns: "urn:zimbra",
			DistributionListActionRequest: reqs
		}
	};
	var respCallback = this._modifyDlResponseHandler.bind(this, fileAsChanged || emailChanged); //there's some issue with fileAsChanged so adding the emailChanged to be on safe side
	appCtxt.getAppController().sendRequest({jsonObj: jsonObj, asyncMode: true, callback: respCallback});

};

ZmContact.prototype._getAddOrRemoveReq =
function(members, add) {
	var req = {
		_jsns: "urn:zimbraAccount",
		dl: {by: "name",
			 _content: this.getEmail()
		},
		action: {
			op: add ? "addMembers" : "removeMembers",
			dlm: []
		}
	};
	for (var i = 0; i < members.length; i++) {
		var member = members[i];
		req.action.dlm.push({_content: member.email});
	}
	return req;

};


ZmContact.prototype._getRenameDlReq =
function(name) {
	return {
		_jsns: "urn:zimbraAccount",
		dl: {by: "name",
			 _content: this.getEmail()
		},
		action: {
			op: "rename",
			newName: {_content: name}
		}
	};
};

ZmContact.prototype._getSetOwnersReq =
function(owners) {
	var ownersPart = [];
	for (var i = 0; i < owners.length; i++) {
		ownersPart.push({
			type: ZmGroupView.GRANTEE_TYPE_USER,
			by: "name",
			_content: owners[i]
		});
	}
	return {
		_jsns: "urn:zimbraAccount",
		dl: {by: "name",
			 _content: this.getEmail()
		},
		action: {
			op: "setOwners",
			owner: ownersPart
		}
	};
};

ZmContact.prototype._getSetMailPolicyReq =
function(mailPolicy, specificMailers) {
	var grantees = [];
	if (mailPolicy == ZmGroupView.MAIL_POLICY_SPECIFIC) {
		for (var i = 0; i < specificMailers.length; i++) {
			grantees.push({
				type: ZmGroupView.GRANTEE_TYPE_EMAIL,
				by: "name",
				_content: specificMailers[i]
			});
		}
	}
	else if (mailPolicy == ZmGroupView.MAIL_POLICY_ANYONE) {
		grantees.push({
			type: ZmGroupView.GRANTEE_TYPE_PUBLIC
		});
	}
	else if (mailPolicy == ZmGroupView.MAIL_POLICY_INTERNAL) {
		grantees.push({
			type: ZmGroupView.GRANTEE_TYPE_ALL
		});
	}
	else if (mailPolicy == ZmGroupView.MAIL_POLICY_MEMBERS) {
		grantees.push({
			type: ZmGroupView.GRANTEE_TYPE_GROUP,
			by: "name",
			_content: this.getEmail()
		});
	}
	else {
		throw "invalid mailPolicy value " + mailPolicy;
	}

	return {
		_jsns: "urn:zimbraAccount",
		dl: {by: "name",
			 _content: this.getEmail()
		},
		action: {
			op: "setRights",
			right: {
				right: "sendToDistList",
				grantee: grantees
			}
		}
	};

};

ZmContact.prototype._addDlAttribute =
function(attrs, mods, name, soapAttrName) {
	var attr = mods[name];
	if (attr === undefined) {
		return;
	}
	attrs.push({n: soapAttrName, _content: attr});
};

ZmContact.prototype._getDlAttributes =
function(mods) {
	var attrs = [];
	this._addDlAttribute(attrs, mods, ZmContact.F_dlDisplayName, "displayName");
	this._addDlAttribute(attrs, mods, ZmContact.F_dlDesc, "description");
	this._addDlAttribute(attrs, mods, ZmContact.F_dlNotes, "zimbraNotes");
	this._addDlAttribute(attrs, mods, ZmContact.F_dlHideInGal, "zimbraHideInGal");
	this._addDlAttribute(attrs, mods, ZmContact.F_dlSubscriptionPolicy, "zimbraDistributionListSubscriptionPolicy");
	this._addDlAttribute(attrs, mods, ZmContact.F_dlUnsubscriptionPolicy, "zimbraDistributionListUnsubscriptionPolicy");

	return attrs;
};


ZmContact.prototype._getCreateDlReq =
function(attr) {
	return {
		_jsns: "urn:zimbraAccount",
		name: attr[ZmContact.F_email],
		a: this._getDlAttributes(attr),
		dynamic: false
	};
};

ZmContact.prototype._getModifyDlAttributesReq =
function(attr) {
	var modAttrs = this._getDlAttributes(attr);
	if (modAttrs.length == 0) {
		return null;
	}
	return {
		_jsns: "urn:zimbraAccount",
		dl: {by: "name",
			 _content: this.getEmail()
		},
		action: {
			op: "modify",
			a: modAttrs
		}
	};
};

ZmContact.prototype._modifyDlResponseHandler =
function(fileAsChanged, result) {
	if (this._handleErrorDl(result)) {
		return;
	}
	appCtxt.setStatusMsg(ZmMsg.dlSaved);

	//for DLs we reload from the server since the server does not send notifications.
	this.clearDlInfo();

	var details = {
		fileAsChanged: fileAsChanged
	};

	this._popView(fileAsChanged);

	this._notify(ZmEvent.E_MODIFY, details);
};

ZmContact.prototype._createDlResponseHandler =
function(result) {
	if (this._handleErrorDl(result, true)) {
		this.attr = {}; //since above in _createDl, we set it to new values prematurely. which would affect next gathering of modified attributes.
		return;
	}
	appCtxt.setStatusMsg(ZmMsg.distributionListCreated);

	this._popView(true);
};

ZmContact.prototype._popView =
function(updateDlList) {
	var controller = AjxDispatcher.run("GetContactController");
	controller.popView(true);
	if (!updateDlList) {
		return;
	}
	var clc = AjxDispatcher.run("GetContactListController");
	if (clc.getFolderId() != ZmFolder.ID_DLS) {
		return;
	}
	ZmAddrBookTreeController.dlFolderClicked(); //This is important in case of new DL created OR a renamed DL, so it would reflect in the list.
};

ZmContact.prototype._handleErrorDl =
function(result, creation) {
	if (!result) {
		return false;
	}
	var batchResp = result.getResponse().BatchResponse;
	var faults = batchResp.Fault;
	if (!faults) {
		return false;
	}
	var ex = ZmCsfeCommand.faultToEx(faults[0]);
	var controller = AjxDispatcher.run("GetContactController");
	controller.popupErrorDialog(creation ? ZmMsg.dlCreateFailed : ZmMsg.dlModifyFailed, ex);
	return true;

};

ZmContact.prototype.clearDlInfo =
function () {
	this.dlMembers = null;
	this.dlInfo = null;
	var app = appCtxt.getApp(ZmApp.CONTACTS);
	app.cacheDL(this.getEmail(), null); //clear the cache for this DL.
	appCtxt.cacheRemove(this.getId()); //also some other cache.
};

/**
 * @private
 */
ZmContact.prototype._handleResponseModify =
function(attr, callback, isAutoSave, result) {
	var resp = result.getResponse().ModifyContactResponse;
	var cn = resp ? resp.cn[0] : null;
	var id = cn ? cn.id : null;
	var groupMembers = cn ? cn.m : null;
	if (groupMembers) {
		this.attr[ZmContact.F_groups] = groupMembers;
		cn._attrs[ZmContact.F_groups] = groupMembers;	
	}

	if (id && id == this.id) {
		if (!isAutoSave) {
			appCtxt.setStatusMsg(this.isGroup() ? ZmMsg.groupSaved : ZmMsg.contactSaved);
		}
		// was this contact moved to another folder?
		if (attr[ZmContact.F_folderId] && this.folderId != attr[ZmContact.F_folderId]) {
			this._setFolder(attr[ZmContact.F_folderId]);
		}
		appCtxt.getApp(ZmApp.CONTACTS).updateIdHash(cn, false);
	} else {
        var detail = ZmMsg.errorTryAgain + "\n" + ZmMsg.errorContact;
        appCtxt.getAppController().setStatusMsg(ZmMsg.errorModifyContact, ZmStatusView.LEVEL_CRITICAL, detail);
	}
	// NOTE: we no longer process callbacks here since notification handling
	//       takes care of everything
};

/**
 * @private
 */
ZmContact.prototype._handleResponseMove =
function(newFolderId, resp) {
	var newFolder = newFolderId && appCtxt.getById(newFolderId);
	var count = 1;
	if (newFolder) {
		appCtxt.setStatusMsg(ZmList.getActionSummary({
			actionTextKey:  'actionMove',
			numItems:       count,
			type:           ZmItem.CONTACT,
			actionArg:      newFolder.name
		}));
	}

	this._notify(ZmEvent.E_MODIFY, resp);
};

/**
 * @private
 */
ZmContact.prototype._setFolder =
function(newFolderId) {
	var folder = appCtxt.getById(this.folderId);
	var fId = folder ? folder.nId : null;
	if (fId == newFolderId) { return; }

	// moving out of a share or into one is handled differently (create then hard delete)
	var newFolder = appCtxt.getById(newFolderId);
	if (this.isShared() || (newFolder && newFolder.link)) {
		if (this.list) {
			this.list.moveItems({items:[this], folder:newFolder});
		}
	} else {
		var jsonObj = {ContactActionRequest:{_jsns:"urn:zimbraMail"}};
		jsonObj.ContactActionRequest.action = {id:this.id, op:"move", l:newFolderId};
		var respCallback = new AjxCallback(this, this._handleResponseMove, [newFolderId]);
		var accountName = appCtxt.multiAccounts && appCtxt.accountList.mainAccount.name;
		appCtxt.getAppController().sendRequest({jsonObj:jsonObj, asyncMode:true, callback:respCallback, accountName:accountName});
	}
};

/**
 * @private
 */
ZmContact.prototype.notifyModify =
function(obj, batchMode) {

	var result = ZmItem.prototype.notifyModify.apply(this, arguments);

	var context = window.parentAppCtxt || window.appCtxt;
	context.clearAutocompleteCache(ZmAutocomplete.AC_TYPE_CONTACT);

	if (result) {
		return result;
	}

	// cache old fileAs/fullName before resetting them
	var oldFileAs = this.getFileAs();
	var oldFullName = this.getFullName();
	this._resetCachedFields();

	var oldAttrCache = {};
	if (obj._attrs) {
		// remove attrs that were not returned back from the server
		var oldAttrs = this.getAttrs();
		for (var a in oldAttrs) {
			oldAttrCache[a] = oldAttrs[a];
			if (obj._attrs[a] == null)
				this.removeAttr(a);
		}

		// set attrs returned by server
		for (var a in obj._attrs) {
			this.setAttr(a, obj._attrs[a]);
		}
		if (obj.m) {
			this.setAttr(ZmContact.F_groups, obj.m);
		}
	}

	var details = {
		attr: obj._attrs,
		oldAttr: oldAttrCache,
		fullNameChanged: (this.getFullName() != oldFullName),
		fileAsChanged: (this.getFileAs() != oldFileAs),
		contact: this
	};

	// update this contact's list per old/new attrs
	for (var listId in this._list) {
		var list = listId && appCtxt.getById(listId);
		if (!list) { continue; }
		list.modifyLocal(obj, details);
	}

	this._notify(ZmEvent.E_MODIFY, obj);
};

/**
 * @private
 */
ZmContact.prototype.notifyDelete =
function() {
	ZmItem.prototype.notifyDelete.call(this);
	var context = window.parentAppCtxt || window.appCtxt;
	context.clearAutocompleteCache(ZmAutocomplete.AC_TYPE_CONTACT);
};

/**
 * Initializes this contact using an email address.
 *
 * @param {AjxEmailAddress|String}	email		an email address or an email string
 * @param {Boolean}	strictName	if <code>true</code>, do not try to set name from user portion of address
 */
ZmContact.prototype.initFromEmail =
function(email, strictName) {
	if (email instanceof AjxEmailAddress) {
		this.setAttr(ZmContact.F_email, email.getAddress());
		this._initFullName(email, strictName);
	} else {
		this.setAttr(ZmContact.F_email, email);
	}
};

/**
 * Initializes this contact using a phone number.
 *
 * @param {String}	phone		the phone string
 * @param {String}	field		the field or company phone if <code>null</code>
 */
ZmContact.prototype.initFromPhone =
function(phone, field) {
	this.setAttr(field || ZmContact.F_companyPhone, phone);
};

/**
 * Gets the email address.
 * 
 * @param {boolean}		asObj	if true, return an AjxEmailAddress
 * 
 * @return	the email address
 */
ZmContact.prototype.getEmail =
function(asObj) {

	var email = (this.getAttr(ZmContact.F_email) ||
				 this.getAttr(ZmContact.F_workEmail1) ||
				 this.getAttr(ZmContact.F_email2) ||
				 this.getAttr(ZmContact.F_workEmail2) ||
				 this.getAttr(ZmContact.F_email3) ||
				 this.getAttr(ZmContact.F_workEmail3));
	
	if (asObj) {
		email = AjxEmailAddress.parse(email);
        if(email){
		    email.isGroup = this.isGroup();
		    email.canExpand = this.canExpand;
        }
	}
	
	return email;
};

/**
 * Returns user's phone number
 * @return {String} phone number
 */
ZmContact.prototype.getPhone = 
function() {
	var phone = (this.getAttr(ZmContact.F_mobilePhone) ||
				this.getAttr(ZmContact.F_workPhone) || 
				this.getAttr(ZmContact.F_homePhone) ||
				this.getAttr(ZmContact.F_otherPhone));
	return phone;
};

    
/**
 * Gets the lookup email address, when an contact object is located using email address we store
 * the referred email address in this variable for easy lookup
 *
 * @param {boolean}		asObj	if true, return an AjxEmailAddress
 *
 * @return	the lookup address
 */
ZmContact.prototype.getLookupEmail =
function(asObj) {
    var email = this._lookupEmail;

    if (asObj && email) {
        email = AjxEmailAddress.parse(email);
        email.isGroup = this.isGroup();
        email.canExpand = this.canExpand;
    }

	return  email;
};

/**
 * Gets the emails.
 * 
 * @return	{Array}	 an array of all valid emails for this contact
 */
ZmContact.prototype.getEmails =
function() {
	var emails = [];
	var attrs = this.getAttrs();
	for (var index = 0; index < ZmContact.EMAIL_FIELDS.length; index++) {
		var field = ZmContact.EMAIL_FIELDS[index];
		for (var i = 1; true; i++) {
			var aname = ZmContact.getAttributeName(field, i);
			if (!attrs[aname]) break;
			emails.push(attrs[aname]);
		}
	}
	return emails;
};

/**
 * Gets the full name.
 * 
 * @return	{String}	the full name
 */
ZmContact.prototype.getFullName =
function(html) {
    var fullNameHtml = null;
	if (!this._fullName || html) {
		var fullName = this.getAttr(ZmContact.X_fullName); // present if GAL contact
		if (fullName) {
			this._fullName = (fullName instanceof Array) ? fullName[0] : fullName;
		}
        else {
            this._fullName = this.getFullNameForDisplay(false);
        }

        if (html) {
            fullNameHtml = this.getFullNameForDisplay(html);
        }
	}

	// as a last resort, set it to fileAs
	if (!this._fullName) {
		this._fullName = this.getFileAs();
	}

	return fullNameHtml || this._fullName;
};

/*
* Gets the fullname for display -- includes (if applicable): prefix, first, middle, maiden, last, suffix
*
* @param {boolean}  if phonetic fields should be used
* @return {String}  the fullname for display
*/
ZmContact.prototype.getFullNameForDisplay =
function(html){
	if (this.isDistributionList()) {
		//I'm not sure where that fullName is set sometime to the display name. This is so complicated
		// I'm trying to set attr[ZmContact.F_dlDisplayName] to the display name but in soem cases it's not.
		return this.getAttr(ZmContact.F_dlDisplayName) || this.getAttr("fullName");
	}
    var prefix = this.getAttr(ZmContact.F_namePrefix);
    var first = this.getAttr(ZmContact.F_firstName);
    var middle = this.getAttr(ZmContact.F_middleName);
    var maiden = this.getAttr(ZmContact.F_maidenName);
    var last = this.getAttr(ZmContact.F_lastName);
    var suffix = this.getAttr(ZmContact.F_nameSuffix);
    var pattern = ZmMsg.fullname;
    if (suffix) {
        pattern = maiden ? ZmMsg.fullnameMaidenSuffix : ZmMsg.fullnameSuffix;
    }
    else if (maiden) {
        pattern = ZmMsg.fullnameMaiden;
    }
    if (appCtxt.get(ZmSetting.LOCALE_NAME) === "ja") {
        var fileAsId = this.getAttr(ZmContact.F_fileAs);
        if (!AjxUtil.isEmpty(fileAsId) && fileAsId !== "1" && fileAsId !== "4" && fileAsId !== "6") {
            /* When Japanese locale is selected, in the most every case, the name should be
             * displayed as "Last First" which is set by the default pattern (ZmMsg_ja.fullname).
             * But if the contact entry's fileAs field explicitly specifies the display
             * format as "First Last", we should override the pattern to lay it out so.
             * For other locales, it is not necessary to override the pattern: The default pattern is
             * already set as "First Last", and even the FileAs specifies as "Last, First", the display
             * name is always expected to be displayed as "First Last".
             */
            pattern = "{0} {1} {2} {4}";
        }
    }
    var formatter = new AjxMessageFormat(pattern);
    var args = [prefix,first,middle,maiden,last,suffix];
    if (!html){
        return AjxStringUtil.trim(formatter.format(args), true);
    }

    return this._getFullNameHtml(formatter, args);
};

/**
 * @param formatter
 * @param parts {Array} Name parts: [prefix,first,middle,maiden,last,suffix]
 */
ZmContact.prototype._getFullNameHtml = function(formatter, parts) {
    var a = [];
    var segments = formatter.getSegments();
    for (var i = 0; i < segments.length; i++) {
        var segment = segments[i];
        if (segment instanceof AjxFormat.TextSegment) {
            a.push(segment.format());
            continue;
        }
        // NOTE: Assume that it's a AjxMessageFormat.MessageSegment
        // NOTE: if not a AjxFormat.TextSegment.
        var index = segment.getIndex();
        var base = parts[index];
        var text = ZmContact.__RUBY_FIELDS[index] && this.getAttr(ZmContact.__RUBY_FIELDS[index]);
        a.push(AjxStringUtil.htmlRubyEncode(base, text));
    }
    return a.join("");
};
ZmContact.__RUBY_FIELDS = [
    null, ZmContact.F_phoneticFirstName, null, null,
    ZmContact.F_phoneticLastName, null
];

/**
 * Gets the tool tip for this contact.
 * 
 * @param	{String}	email		the email address
 * @param	{Boolean}	isGal		(not used)
 * @param	{String}	hint		the hint text
 * @return	{String}	the tool tip in HTML
 */
ZmContact.prototype.getToolTip =
function(email, isGal, hint) {
	// XXX: we dont cache tooltip info anymore since its too dynamic :/
	// i.e. IM status can change anytime so always rebuild tooltip and bug 13834
	var subs = {
		contact: this,
		entryTitle: this.getFileAs(),
		hint: hint
	};

	return (AjxTemplate.expand("abook.Contacts#Tooltip", subs));
};

/**
 * Gets the filing string for this contact, computing it if necessary.
 * 
 * @param	{Boolean}	lower		<code>true</code> to use lower case
 * @return	{String}	the file as string
 */
ZmContact.prototype.getFileAs =
function(lower) {
	// update/null if modified
	if (!this._fileAs) {
		this._fileAs = ZmContact.computeFileAs(this);
		this._fileAsLC = this._fileAs ? this._fileAs.toLowerCase() : null;
	}
	// if for some reason fileAsLC is not set even though fileAs is, reset it
	if (lower && !this._fileAsLC) {
		this._fileAsLC = this._fileAs.toLowerCase();
	}
	return lower ? this._fileAsLC : this._fileAs;
};

/**
 * Gets the filing string for this contact, from the email address (used in case no name exists).
 * todo - maybe return this from getFileAs, but there are a lot of callers to getFileAs, and not sure
 * of the implications on all the use-cases.
 *
 * @return	{String}	the file as string
 */
ZmContact.prototype.getFileAsNoName = function() {
	return [ZmMsg.noName, this.getEmail()].join(" ");
};

/**
 * Gets the header.
 * 
 * @return	{String}	the header
 */
ZmContact.prototype.getHeader =
function() {
	return this.id ? this.getFileAs() : ZmMsg.newContact;
};

ZmContact.NO_MAX_IMAGE_WIDTH = ZmContact.NO_MAX_IMAGE_HEIGHT = - 1;

/**
 * Get the image URL.
 *
 * Please note that maxWidth and maxHeight are hints, as they have no
 * effect on Zimlet-supplied images.
 *
 * maxWidth {int} max pixel width (optional - default 48, or pass ZmContact.NO_MAX_IMAGE_WIDTH if full size image is required)
 * maxHeight {int} max pixel height (optional - default to maxWidth, or pass ZmContact.NO_MAX_IMAGE_HEIGHT if full size image is required)
 * @return	{String}	the image URL
 */
ZmContact.prototype.getImageUrl =
function(maxWidth, maxHeight) {
  	var image = this.getAttr(ZmContact.F_image);
	var imagePart  = image && image.part || this.getAttr(ZmContact.F_imagepart); //see bug 73146

	if (!imagePart) {
		return this.getAttr(ZmContact.F_zimletImage);  //return zimlet populated image only if user-uploaded image is not there.
	}
  	var msgFetchUrl = appCtxt.get(ZmSetting.CSFE_MSG_FETCHER_URI);
	var maxWidthStyle = "";
	if (maxWidth !== ZmContact.NO_MAX_IMAGE_WIDTH) {
		maxWidth = maxWidth || 48;
		maxWidthStyle = ["&max_width=", maxWidth].join("");
	}
	var maxHeightStyle = "";
	if (maxHeight !== ZmContact.NO_MAX_IMAGE_HEIGHT) {
		maxHeight = maxHeight ||
			(maxWidth !== ZmContact.NO_MAX_IMAGE_WIDTH ? maxWidth : 48);
		maxHeightStyle = ["&max_height=", maxHeight].join("");
	}
  	return  [msgFetchUrl, "&id=", this.id, "&part=", imagePart, maxWidthStyle, maxHeightStyle, "&t=", (new Date()).getTime()].join("");
};

ZmContact.prototype.addModifyZimletImageToBatch =
function(batchCmd, image) {
	var attr = {};
	if (this.getAttr(ZmContact.F_zimletImage) === image) {
		return; //no need to update if same
	}
	attr[ZmContact.F_zimletImage] = image;
	batchCmd.add(this.modify.bind(this, attr, null, true));
};

/**
 * Gets the company field. Company field has a getter b/c fileAs may be the Company name so
 * company field should return "last, first" name instead *or* prepend the title
 * if fileAs is not Company (assuming it exists).
 * 
 * @return	{String}	the company
 */
ZmContact.prototype.getCompanyField =
function() {

	var attrs = this.getAttrs();
	if (attrs == null) return null;

	var fa = parseInt(attrs.fileAs);
	var val = [];
	var idx = 0;

	if (fa == ZmContact.FA_LAST_C_FIRST || fa == ZmContact.FA_FIRST_LAST) {
		// return the title, company name
		if (attrs.jobTitle) {
			val[idx++] = attrs.jobTitle;
			if (attrs.company)
				val[idx++] = ", ";
		}
		if (attrs.company)
			val[idx++] = attrs.company;

	} else if (fa == ZmContact.FA_COMPANY) {
		// return the first/last name
		if (attrs.lastName) {
			val[idx++] = attrs.lastName;
			if (attrs.firstName)
				val[idx++] = ", ";
		}

		if (attrs.firstName)
			val[idx++] = attrs.firstName;

		if (attrs.jobTitle)
			val[idx++] = " (" + attrs.jobTitle + ")";

	} else {
		// just return the title
		if (attrs.jobTitle) {
			val[idx++] = attrs.jobTitle;
			// and/or company name if applicable
			if (attrs.company && (attrs.fileAs == null || fa == ZmContact.FA_LAST_C_FIRST || fa == ZmContact.FA_FIRST_LAST))
				val[idx++] = ", ";
		}
		if (attrs.company && (attrs.fileAs == null || fa == ZmContact.FA_LAST_C_FIRST || fa == ZmContact.FA_FIRST_LAST))
			 val[idx++] = attrs.company;
	}
	if (val.length == 0) return null;
	return val.join("");
};

/**
 * Gets the work address.
 * 
 * @param	{Object}	instance		(not used)
 * @return	{String}	the work address
 */
ZmContact.prototype.getWorkAddrField =
function(instance) {
	var attrs = this.getAttrs();
	return this._getAddressField(attrs.workStreet, attrs.workCity, attrs.workState, attrs.workPostalCode, attrs.workCountry);
};

/**
 * Gets the home address.
 * 
 * @param	{Object}	instance		(not used)
 * @return	{String}	the home address
 */
ZmContact.prototype.getHomeAddrField =
function(instance) {
	var attrs = this.getAttrs();
	return this._getAddressField(attrs.homeStreet, attrs.homeCity, attrs.homeState, attrs.homePostalCode, attrs.homeCountry);
};

/**
 * Gets the other address.
 * 
 * @param	{Object}	instance		(not used)
 * @return	{String}	the other address
 */
ZmContact.prototype.getOtherAddrField =
function(instance) {
	var attrs = this.getAttrs();
	return this._getAddressField(attrs.otherStreet, attrs.otherCity, attrs.otherState, attrs.otherPostalCode, attrs.otherCountry);
};

/**
 * Gets the address book.
 * 
 * @return	{ZmAddrBook}	the address book
 */
ZmContact.prototype.getAddressBook =
function() {
	if (!this.addrbook) {
		this.addrbook = appCtxt.getById(this.folderId);
	}
	return this.addrbook;
};

/**
 * @private
 */
ZmContact.prototype._getAddressField =
function(street, city, state, zipcode, country) {
	if (street == null && city == null && state == null && zipcode == null && country == null) return null;

	var html = [];
	var idx = 0;

	if (street) {
		html[idx++] = street;
		if (city || state || zipcode)
			html[idx++] = "\n";
	}

	if (city) {
		html[idx++] = city;
		if (state)
			html[idx++] = ", ";
		else if (zipcode)
			html[idx++] = " ";
	}

	if (state) {
		html[idx++] = state;
		if (zipcode)
			html[idx++] = " ";
	}

	if (zipcode)
		html[idx++] = zipcode;

	if (country)
		html[idx++] = "\n" + country;

	return html.join("");
};

/**
 * Sets the full name based on an email address.
 * 
 * @private
 */
ZmContact.prototype._initFullName =
function(email, strictName) {
	var name = email.getName();
	name = AjxStringUtil.trim(name.replace(AjxEmailAddress.commentPat, '')); // strip comment (text in parens)

	if (name && name.length) {
		this._setFullName(name, [" "]);
	} else if (!strictName) {
		name = email.getAddress();
		if (name && name.length) {
			var i = name.indexOf("@");
			if (i == -1) return;
			name = name.substr(0, i);
			this._setFullName(name, [".", "_"]);
		}
	}
};

/**
 * Tries to extract a set of name components from the given text, with the
 * given list of possible delimiters. The first delimiter contained in the
 * text will be used. If none are found, the first delimiter in the list is used.
 * 
 * @private
 */
ZmContact.prototype._setFullName =
function(text, delims) {
	var delim = delims[0];
	for (var i = 0; i < delims.length; i++) {
		if (text.indexOf(delims[i]) != -1) {
			delim = delims[i];
			break;
		}
	}
    var parts = text.split(delim);
    var func = this["__setFullName_"+AjxEnv.DEFAULT_LOCALE] || this.__setFullName;
    func.call(this, parts, text, delims);
};

ZmContact.prototype.__setFullName = function(parts, text, delims) {
    this.setAttr(ZmContact.F_firstName, parts[0]);
    if (parts.length == 2) {
        this.setAttr(ZmContact.F_lastName, parts[1]);
    } else if (parts.length == 3) {
        this.setAttr(ZmContact.F_middleName, parts[1]);
        this.setAttr(ZmContact.F_lastName, parts[2]);
    }
};
ZmContact.prototype.__setFullName_ja = function(parts, text, delims) {
    if (parts.length > 2) {
        this.__setFullName(parts, text, delims);
        return;
    }
    // TODO: Perhaps do some analysis to auto-detect Japanese vs.
    // TODO: non-Japanese names. For example, if the name text is
    // TODO: comprised of kanji, treat it as "last first"; else if
    // TODO: first part is all uppercase, treat it as "last first";
    // TODO: else treat it as "first last".
    this.setAttr(ZmContact.F_lastName, parts[0]);
    if (parts.length > 1) {
        this.setAttr(ZmContact.F_firstName, parts[1]);
    }
};
ZmContact.prototype.__setFullName_ja_JP = ZmContact.prototype.__setFullName_ja;

/**
 * @private
 */
ZmContact.prototype._addRequestAttr =
function(cn, name, value) {
	var a = {n:name};
	if (name == ZmContact.F_image && AjxUtil.isString(value) && value.length) {
		// handle contact photo
		if (value.indexOf("aid_") != -1) {
			a.aid = value.substring(4);
		} else {
			a.part = value.substring(5);
		}
	} else {
		a._content = value || "";
	}

    if (value instanceof Array) {
        if (!cn._attrs)
            cn._attrs = {};
        cn._attrs[name] = value || "";
    }
    else  {
        if (!cn.a)
            cn.a = [];
        cn.a.push(a);
    }
};
	
ZmContact.prototype._addContactGroupAttr = 
function(cn, group) {
	var groupMembers = group[ZmContact.F_groups];
	for (var i = 0; i < groupMembers.length; i++) {
		var member = groupMembers[i];
		if (!cn.m) {
			cn.m = [];
		}

		var m = {type: member.type,	value: member.value}; //for the JSON object this is all we need.
		if (member.op) {
			m.op = member.op; //this is only for modify, not for create.
		}
		cn.m.push(m);
	}
};

/**
 * Reset computed fields.
 * 
 * @private
 */
ZmContact.prototype._resetCachedFields =
function() {
	this._fileAs = this._fileAsLC = this._fullName = null;
};

/**
 * Parse contact node.
 * 
 * @private
 */
ZmContact.prototype._loadFromDom =
function(node) {
	this.isLoaded = true;
	this.rev = node.rev;
	this.sf = node.sf || node._attrs.sf;
	if (!this.isGal) {
		this.folderId = node.l;
	}
	this.created = node.cd;
	this.modified = node.md;

	this.attr = node._attrs || {};
	if (node.m) {
		this.attr[ZmContact.F_groups] = node.m;
	}

	this.ref = node.ref || this.attr.dn; //bug 78425
	
	// for shared contacts, we get these fields outside of the attr part
	if (node.email)		{ this.attr[ZmContact.F_email] = node.email; }
	if (node.email2)	{ this.attr[ZmContact.F_email2] = node.email2; }
	if (node.email3)	{ this.attr[ZmContact.F_email3] = node.email3; }

	// in case attrs are coming in from an external GAL, make an effort to map them, including multivalued ones
	this.attr = ZmContact.mapAttrs(this.attr);

    //the attr groups is returned as [] so check both null and empty array to set the type
    var groups = this.attr[ZmContact.F_groups];
    if(!groups || (groups instanceof Array && groups.length == 0)) {
        this.type = ZmItem.CONTACT;
    }
    else {
        this.type = ZmItem.GROUP;
    }

	// check if the folderId is found in our address book (otherwise, we assume
	// this contact to be a shared contact)
	var ac = window.parentAppCtxt || window.appCtxt;
	this.addrbook = ac.getById(this.folderId);

	this._parseTagNames(node.tn);

	// dont process flags for shared contacts until we get server support
	if (!this.isShared()) {
		this._parseFlags(node.f);
	} else {
		// shared contacts are never fully loaded since we never cache them
		this.isLoaded = false;
	}

	// bug: 22174
	// We ignore the server's computed file-as property and instead
	// format it based on the user's locale.
	this._fileAs = ZmContact.computeFileAs(this);

	// Is this a distribution list?
	this.isDL = this.isDistributionList();
	if (this.isDL) {
		this.dlInfo = { //this is minimal DL info, available mainly to allow to know whether to show the lock or not.
			isMinimal: true,
			isMember: node.isMember,
			isOwner: node.isOwner,
			subscriptionPolicy: this.attr.zimbraDistributionListSubscriptionPolicy,
			unsubscriptionPolicy: this.attr.zimbraDistributionListUnsubscriptionPolicy,
			displayName: node.d || "",
			hideInGal: this.attr.zimbraHideInGal == "TRUE"
		};

		this.canExpand = node.exp !== false; //default to true, since most cases this is implicitly true if not returned. See bug 94867
		var emails = this.getEmails();
		var ac = window.parentAppCtxt || window.appCtxt;
		for (var i = 0; i < emails.length; i++) {
			ac.setIsExpandableDL(emails[i], this.canExpand);
		}
	}
	//S/MIME: If user certificate is present, include it into contact's object
	if (node.certificate) {
		this.certificate = node.certificate;
	}
};

/**
 * Gets display text for an attendee. Prefers name over email.
 *
 * @param {constant}	type		the attendee type
 * @param {Boolean}	shortForm		if <code>true</code>, return only name or email
 * @return	{String}	the attendee
 */
ZmContact.prototype.getAttendeeText =
function(type, shortForm) {
	var email = this.getEmail(true);
	return (email?email.toString(shortForm || (type && type != ZmCalBaseItem.PERSON)):"");
};

/**
 * Gets display text for an attendee. Prefers name over email.
 *
 * @param {constant}	type		the attendee type
 * @param {Boolean}	shortForm		if <code>true</code>, return only name or email
 * @return	{String}	the attendee
 */
ZmContact.prototype.getAttendeeKey =
function() {
	var email = this.getLookupEmail() || this.getEmail();
	var name = this.getFullName();
	return email ? email : name;
};

/**
 * Gets the unknown fields.
 * 
 * @param	{function}	[sortByNameFunc]	sort by function
 * @return	{Array}	an array of field name/value pairs
 */
ZmContact.prototype.getUnknownFields = function(sortByNameFunc) {
	var map = ZmContact.__FIELD_MAP;
	if (!map) {
		map = ZmContact.__FIELD_MAP = {};
		for (var i = 0; i < ZmContact.DISPLAY_FIELDS; i++) {
			map[ZmContact.DISPLAY_FIELDS[i]] = true;
		}
	}
	var fields = [];
	var attrs = this.getAttrs();
	for (var aname in attrs) {
		var field = ZmContact.getPrefix(aname);
		if (map[aname]) continue;
		fields.push(field);
	}
	return this.getFields(fields, sortByNameFunc);
};

/**
 * Gets the fields.
 * 
 * @param	{Array}	field		the fields
 * @param	{function}	[sortByNameFunc]	sort by function
 * @return	{Array}	an array of field name/value pairs
 */
ZmContact.prototype.getFields =
function(fields, sortByNameFunc) {
	// TODO: [Q] Should sort function handle just the field names or the attribute names?
	var selection;
	var attrs = this.getAttrs();
	for (var index = 0; index < fields.length; index++) {
		for (var i = 1; true; i++) {
			var aname = ZmContact.getAttributeName(fields[index], i);
			if (!attrs[aname]) break;
			if (!selection) selection = {};
			selection[aname] = attrs[aname];
		}
	}
	if (sortByNameFunc && selection) {
		var keys = AjxUtil.keys(selection);
		keys.sort(sortByNameFunc);
		var nfields = {};
		for (var i = 0; i < keys; i++) {
			var key = keys[i];
			nfields[key] = fields[key];
		}
		selection = nfields;
	}
	return selection;
};

/**
 * Returns a list of distribution list members for this contact. Only the
 * requested range is returned.
 *
 * @param offset	{int}			offset into list to start at
 * @param limit		{int}			number of members to fetch and return
 * @param callback	{AjxCallback}	callback to run with results
 */
ZmContact.prototype.getDLMembers =
function(offset, limit, callback) {

	var result = {list:[], more:false, isDL:{}};
	if (!this.isDL) { return result; }

	var email = this.getEmail();
	var app = appCtxt.getApp(ZmApp.CONTACTS);
	var dl = app.getDL(email);
	if (!dl) {
		dl = result;
		dl.more = true;
		app.cacheDL(email, dl);
	}

	limit = limit || ZmContact.DL_PAGE_SIZE;
	var start = offset || 0;
	var end = (offset + limit) - 1;

	// see if we already have the requested members, or know that we don't
	if (dl.list.length >= end + 1 || !dl.more) {
		var list = dl.list.slice(offset, end + 1);
		result = {list:list, more:dl.more || (dl.list.length > end + 1), isDL:dl.isDL};
		DBG.println("dl", "found cached DL members");
		this._handleResponseGetDLMembers(start, limit, callback, result);
		return;
	}

	DBG.println("dl", "server call " + offset + " / " + limit);
	if (!dl.total || (offset < dl.total)) {
		var jsonObj = {GetDistributionListMembersRequest:{_jsns:"urn:zimbraAccount", offset:offset, limit:limit}};
		var request = jsonObj.GetDistributionListMembersRequest;
		request.dl = {_content: this.getEmail()};
		var respCallback = new AjxCallback(this, this._handleResponseGetDLMembers, [offset, limit, callback]);
		appCtxt.getAppController().sendRequest({jsonObj:jsonObj, asyncMode:true, callback:respCallback});
	} else {
		this._handleResponseGetDLMembers(start, limit, callback, result);
	}
};

ZmContact.prototype._handleResponseGetDLMembers =
function(offset, limit, callback, result, resp) {

	if (resp || !result.list) {
		var list = [];
		resp = resp || result.getResponse();  //if response is passed, take it. Otherwise get it from result
		resp = resp.GetDistributionListMembersResponse;
		var dl = appCtxt.getApp(ZmApp.CONTACTS).getDL(this.getEmail());
		var more = dl.more = resp.more;
		var isDL = {};
		var members = resp.dlm;
		if (members && members.length) {
			for (var i = 0, len = members.length; i < len; i++) {
				var member = members[i]._content;
				list.push(member);
				dl.list[offset + i] = member;
				if (members[i].isDL) {
					isDL[member] = dl.isDL[member] = true;
				}
			}
		}
		dl.total = resp.total;
		DBG.println("dl", list.join("<br>"));
		var result = {list:list, more:more, isDL:isDL};
	}
	DBG.println("dl", "returning list of " + result.list.length + ", more is " + result.more);
	if (callback) {
		callback.run(result);
	}
	else { //synchronized case - see ZmContact.prototype.getDLMembers above
		return result;
	}
};

/**
 * Returns a list of all the distribution list members for this contact.
 *
 * @param callback	{AjxCallback}	callback to run with results
 */
ZmContact.prototype.getAllDLMembers =
function(callback) {

	var result = {list:[], more:false, isDL:{}};
	if (!this.isDL) { return result; }

	var dl = appCtxt.getApp(ZmApp.CONTACTS).getDL(this.getEmail());
	if (dl && !dl.more) {
		result = {list:dl.list.slice(), more:false, isDL:dl.isDL};
		callback.run(result);
		return;
	}

	var nextCallback = new AjxCallback(this, this._getNextDLChunk, [callback]);
	this.getDLMembers(dl ? dl.list.length : 0, null, nextCallback);
};

ZmContact.prototype._getNextDLChunk =
function(callback, result) {

	var dl = appCtxt.getApp(ZmApp.CONTACTS).getDL(this.getEmail());
	if (result.more) {
		var nextCallback = new AjxCallback(this, this._getNextDLChunk, [callback]);
		this.getDLMembers(dl.list.length, null, nextCallback);
	} else {
		result.list = dl.list.slice();
		callback.run(result);
	}
};

/**
 * Gets the contact from cache handling parsing of contactId
 * 
 * @param contactId {String} contact id
 * @return contact {ZmContact} contact or null
 * @private
 */
ZmContact.getContactFromCache =
function(contactId) {
	var userZid = appCtxt.accountList.mainAccount.id;
	var contact = null;
	if (contactId && contactId.indexOf(userZid + ":") !=-1) {
		//strip off the usersZid to pull from cache
		var arr = contactId.split(userZid + ":");
		contact = arr && arr.length > 1 ? appCtxt.cacheGet(arr[1]) : appCtxt.cacheGet(contactId);
	}
	else {
		contact = appCtxt.cacheGet(contactId);
	}
	if (contact instanceof ZmContact) {
		return contact;
	}
	return null;
};

// For mapAttrs(), prepare a hash where each key is the base name of an attr (without an ending number and lowercased),
// and the value is a numerically sorted list of attr names in their original form.
ZmContact.ATTR_VARIANTS = {};
ZmContact.IGNORE_ATTR_VARIANT = {};
ZmContact.IGNORE_ATTR_VARIANT[ZmContact.F_groups] = true;

ZmContact.initAttrVariants = function(attrClass) {
	var keys = Object.keys(attrClass),
		len = keys.length, key, i, attr,
		attrs = [];

	// first, grab all the attr names
	var ignoreVariant = attrClass.IGNORE_ATTR_VARIANT || {};
	for (i = 0; i < len; i++) {
		key = keys[i];
		if (key.indexOf('F_') === 0) {
			attr = attrClass[key];
			if (!ignoreVariant[attr]) {
				attrs.push(attr);
			}
		}
	}

	// sort numerically, eg so that we get ['email', 'email2', 'email10'] in right order
	var numRegex = /^([a-zA-Z]+)(\d+)$/;
	attrs.sort(function(a, b) {
		var aMatch = a.match(numRegex),
			bMatch = b.match(numRegex);
		// check if both are numbered attrs with same base
		if (aMatch && bMatch && aMatch[1] === bMatch[1]) {
			return aMatch[2] - bMatch[2];
		}
		else {
			return a > b ? 1 : (a < b ? -1 : 0);
		}
	});

	// construct hash mapping generic base name to its iterated attr names
	var attr, base;
	for (i = 0; i < attrs.length; i++) {
		attr = attrs[i];
		base = attr.replace(/\d+$/, '').toLowerCase();
		if (!ZmContact.ATTR_VARIANTS[base]) {
			ZmContact.ATTR_VARIANTS[base] = [];
		}
		ZmContact.ATTR_VARIANTS[base].push(attr);
	}
};
ZmContact.initAttrVariants(ZmContact);

/**
 * Takes a hash of attrs and values and maps it to our attr names as best as it can. Scalar attrs will map if they
 * have the same name or only differ by case. A multivalued attr will map to a set of our attributes that share the
 * same case-insensitive base name. Some examples:
 *
 *      FIRSTNAME: "Mildred"    =>      firstName: "Mildred"
 *      email: ['a', 'b']       =>      email: 'a',
 *                                      email2: 'b'
 *      WorkEmail: ['y', 'z']   =>      workEmail1: 'y',
 *                                      workEmail2: 'z'
 *      IMaddress: ['f', 'g']   =>      imAddress1: 'f',
 *                                      imAddress2: 'g'
 *
 * @param   {Object}    attrs       hash of attr names/values
 *
 * @returns {Object}    hash of attr names/values using known attr names ZmContact.F_*
 */
ZmContact.mapAttrs = function(attrs) {

	var attr, value, baseAttrs, newAttrs = {};
	for (attr in attrs) {
		value = attrs[attr];
		if (value) {
			baseAttrs = ZmContact.ATTR_VARIANTS[attr.toLowerCase()];
			if (baseAttrs) {
				value = AjxUtil.toArray(value);
				var len = Math.min(value.length, baseAttrs.length), i;
				for (i = 0; i < len; i++) {
					newAttrs[baseAttrs[i]] = value[i];
				}
			} else {
				// Any overlooked/ignored attributes are simply passed along
				newAttrs[attr] = value;
			}
		}
	}
	return newAttrs;
};

// these need to be kept in sync with ZmContact.F_*
ZmContact._AB_FIELD = {
	firstName:				ZmMsg.AB_FIELD_firstName,		// file as info
	lastName:				ZmMsg.AB_FIELD_lastName,
	middleName:				ZmMsg.AB_FIELD_middleName,
	fullName:				ZmMsg.AB_FIELD_fullName,
	jobTitle:				ZmMsg.AB_FIELD_jobTitle,
	company:				ZmMsg.AB_FIELD_company,
	department:				ZmMsg.AB_FIELD_department,
	email:					ZmMsg.AB_FIELD_email,			// email addresses
	email2:					ZmMsg.AB_FIELD_email2,
	email3:					ZmMsg.AB_FIELD_email3,
	imAddress1:				ZmMsg.AB_FIELD_imAddress1,		// IM addresses
	imAddress2:				ZmMsg.AB_FIELD_imAddress2,
	imAddress3:				ZmMsg.AB_FIELD_imAddress3,
	image: 					ZmMsg.AB_FIELD_image,			// contact photo
	attachment:				ZmMsg.AB_FIELD_attachment,
	workStreet:				ZmMsg.AB_FIELD_street,			// work address info
	workCity:				ZmMsg.AB_FIELD_city,
	workState:				ZmMsg.AB_FIELD_state,
	workPostalCode:			ZmMsg.AB_FIELD_postalCode,
	workCountry:			ZmMsg.AB_FIELD_country,
	workURL:				ZmMsg.AB_FIELD_URL,
	workPhone:				ZmMsg.AB_FIELD_workPhone,
	workPhone2:				ZmMsg.AB_FIELD_workPhone2,
	workFax:				ZmMsg.AB_FIELD_workFax,
	assistantPhone:			ZmMsg.AB_FIELD_assistantPhone,
	companyPhone:			ZmMsg.AB_FIELD_companyPhone,
	callbackPhone:			ZmMsg.AB_FIELD_callbackPhone,
	homeStreet:				ZmMsg.AB_FIELD_street,			// home address info
	homeCity:				ZmMsg.AB_FIELD_city,
	homeState:				ZmMsg.AB_FIELD_state,
	homePostalCode:			ZmMsg.AB_FIELD_postalCode,
	homeCountry:			ZmMsg.AB_FIELD_country,
	homeURL:				ZmMsg.AB_FIELD_URL,
	homePhone:				ZmMsg.AB_FIELD_homePhone,
	homePhone2:				ZmMsg.AB_FIELD_homePhone2,
	homeFax:				ZmMsg.AB_FIELD_homeFax,
	mobilePhone:			ZmMsg.AB_FIELD_mobilePhone,
	pager:					ZmMsg.AB_FIELD_pager,
	carPhone:				ZmMsg.AB_FIELD_carPhone,
	otherStreet:			ZmMsg.AB_FIELD_street,			// other info
	otherCity:				ZmMsg.AB_FIELD_city,
	otherState:				ZmMsg.AB_FIELD_state,
	otherPostalCode:		ZmMsg.AB_FIELD_postalCode,
	otherCountry:			ZmMsg.AB_FIELD_country,
	otherURL:				ZmMsg.AB_FIELD_URL,
	otherPhone:				ZmMsg.AB_FIELD_otherPhone,
	otherFax:				ZmMsg.AB_FIELD_otherFax,
	notes:					ZmMsg.notes,					// misc fields
	birthday:				ZmMsg.AB_FIELD_birthday
};

ZmContact._AB_FILE_AS = {
	1:						ZmMsg.AB_FILE_AS_lastFirst,
	2:						ZmMsg.AB_FILE_AS_firstLast,
	3:						ZmMsg.AB_FILE_AS_company,
	4:						ZmMsg.AB_FILE_AS_lastFirstCompany,
	5:						ZmMsg.AB_FILE_AS_firstLastCompany,
	6:						ZmMsg.AB_FILE_AS_companyLastFirst,
	7:						ZmMsg.AB_FILE_AS_companyFirstLast
};

} // if (!window.ZmContact)
}
if (AjxPackage.define("zimbraMail.abook.model.ZmContactList")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */

/**
 * @overview
 * This file contains the contact list class.
 * 
 */

/**
 * Create a new, empty contact list.
 * @class
 * This class represents a list of contacts. In general, the list is the result of a
 * search. It may be the result of a <code>&lt;GetContactsRequest&gt;</code>, which returns all of the user's
 * local contacts. That list is considered to be canonical.
 * <p>
 * Loading of all local contacts has been optimized by delaying the creation of {@link ZmContact} objects until
 * they are needed. That has a big impact on IE, and not much on Firefox. Loading a subset
 * of attributes did not have much impact on load time, probably because a large majority
 * of contacts contain only those minimal fields.</p>
 *
 * @author Conrad Damon
 *
 * @param {ZmSearch}	search	the search that generated this list
 * @param {Boolean}	isGal		if <code>true</code>, this is a list of GAL contacts
 * @param {constant}	type		the item type
 * 
 * @extends		ZmList
 */
ZmContactList = function(search, isGal, type) {

	if (arguments.length == 0) { return; }
	type = type || ZmItem.CONTACT;
	ZmList.call(this, type, search);

	this.isGal = (isGal === true);
	this.isCanonical = false;
	this.isLoaded = false;

	this._app = appCtxt.getApp(ZmApp.CONTACTS);
	if (!this._app) { 
		this._emailToContact = this._phoneToContact = {};
		return;
	}
	this._emailToContact = this._app._byEmail;
	this._phoneToContact = this._app._byPhone;

	this._alwaysUpdateHashes = true; // Should we update the phone & IM fast-lookup hashes even when account features don't require it? (bug #60411)
};

ZmContactList.prototype = new ZmList;
ZmContactList.prototype.constructor = ZmContactList;

ZmContactList.prototype.isZmContactList = true;
ZmContactList.prototype.toString = function() { return "ZmContactList"; };




// Constants

// Support for loading user's local contacts from a large string

ZmContactList.URL = "/Contacts";	// REST URL for loading user's local contacts
ZmContactList.URL_ARGS = { fmt: 'cf', t: 2, all: 'all' }; // arguments for the URL above
ZmContactList.CONTACT_SPLIT_CHAR	= '\u001E';	// char for splitting string into contacts
ZmContactList.FIELD_SPLIT_CHAR		= '\u001D';	// char for splitting contact into fields
// fields that belong to a contact rather than its attrs
ZmContactList.IS_CONTACT_FIELD = {"id":true, "l":true, "d":true, "fileAsStr":true, "rev":true};



/**
 * @private
 */
ZmContactList.prototype.addLoadedCallback =
function(callback) {
	if (this.isLoaded) {
		callback.run();
		return;
	}
	if (!this._loadedCallbacks) {
		this._loadedCallbacks = [];
	}
	this._loadedCallbacks.push(callback);
};

/**
 * @private
 */
ZmContactList.prototype._finishLoading =
function() {
	DBG.timePt("done loading " + this.size() + " contacts");
	this.isLoaded = true;
	if (this._loadedCallbacks) {
		var callback;
		while (callback = this._loadedCallbacks.shift()) {
			callback.run();
		}
	}
};

/**
 * Retrieves the contacts from the back end, and parses the response. The list is then sorted.
 * This method is used only by the canonical list of contacts, in order to load their content.
 * <p>
 * Loading a minimal set of attributes did not result in a significant performance gain.
 * </p>
 * 
 * @private
 */
ZmContactList.prototype.load =
function(callback, errorCallback, accountName) {
	// only the canonical list gets loaded
	this.isCanonical = true;
	var respCallback = new AjxCallback(this, this._handleResponseLoad, [callback]);
	DBG.timePt("requesting contact list", true);
    if(appCtxt.isExternalAccount()) {
        //Do not make a call in case of external user
        //The rest url constructed wont exist in case of external user
        if (callback) {
		    callback.run();
	    }
        return;
    }
	var args = ZmContactList.URL_ARGS;

	// bug 74609: suppress overzealous caching by IE
	if (AjxEnv.isIE) {
		args = AjxUtil.hashCopy(args);
		args.sid = ZmCsfeCommand.getSessionId();
	}

	var params = {asyncMode:true, noBusyOverlay:true, callback:respCallback, errorCallback:errorCallback, offlineCallback:callback};
	params.restUri = AjxUtil.formatUrl({
		path:["/home/", (accountName || appCtxt.getUsername()),
	          ZmContactList.URL].join(""),
	    qsArgs: args, qsReset:true
	});
	DBG.println(AjxDebug.DBG1, "loading contacts from " + params.restUri);
	appCtxt.getAppController().sendRequest(params);

	ZmContactList.addDlFolder();
	
};

/**
 * @private
 */
ZmContactList.prototype._handleResponseLoad =
function(callback, result) {
	DBG.timePt("got contact list");
	var text = result.getResponse();
    if (text && typeof text !== 'string'){
        text = text._data;
    }
	var derefList = [];
	if (text) {
		var contacts = text.split(ZmContactList.CONTACT_SPLIT_CHAR);
		var derefBatchCmd = new ZmBatchCommand(true, null, true);
		for (var i = 0, len = contacts.length; i < len; i++) {
			var fields = contacts[i].split(ZmContactList.FIELD_SPLIT_CHAR);
			var contact = {}, attrs = {};
			var groupMembers = [];
			var foundDeref = false;
			for (var j = 0, len1 = fields.length; j < len1; j += 2) {
				if (ZmContactList.IS_CONTACT_FIELD[fields[j]]) {
					contact[fields[j]] = fields[j + 1];
				} else {
					var value = fields[j+1];
					switch (fields[j]) {
						case ZmContact.F_memberC:
							groupMembers.push({type: ZmContact.GROUP_CONTACT_REF, value: value});
							foundDeref = true; //load shared contacts
							break;
						case ZmContact.F_memberG:
							groupMembers.push({type: ZmContact.GROUP_GAL_REF, value: value});
							foundDeref = true;
							break;
						case ZmContact.F_memberI:
							groupMembers.push({type: ZmContact.GROUP_INLINE_REF, value: value});
							foundDeref = true;
							break;
						default:
							attrs[fields[j]] = value;
					}
				}
			}
			if (attrs[ZmContact.F_type] === "group") { //set only for group.
				attrs[ZmContact.F_groups] = groupMembers;
			}
			if (foundDeref) {
				//batch group members for deref loading
				var dummy = new ZmContact(contact["id"], this);
				derefBatchCmd.add(new AjxCallback(dummy, dummy.load, [null, null, derefBatchCmd, true]));
			}
			contact._attrs = attrs;
			this._addContact(contact);
		}
		derefBatchCmd.run();
	}

	this._finishLoading();

	if (callback) {
		callback.run();
	}
};

/**
 * @static
 */
ZmContactList.addDlFolder =
function() {

	if (!appCtxt.get(ZmSetting.DLS_FOLDER_ENABLED)) {
		return;
	}

	var dlsFolder = appCtxt.getById(ZmOrganizer.ID_DLS);

	var root = appCtxt.getById(ZmOrganizer.ID_ROOT);
	if (!root) { return; }

	if (dlsFolder && root.getById(ZmOrganizer.ID_DLS)) {
		//somehow (after a refresh block, can be reprod using $set:refresh. ZmClientCmdHandler.prototype.execute_refresh) the DLs folder object is removed from under the root (but still cached in appCtxt). So making sure it's there.
		return;
	}

	if (!dlsFolder) {
		var params = {
			id: ZmOrganizer.ID_DLS,
			name: ZmMsg.distributionLists,
			parent: root,
			tree: root.tree,
			type: ZmOrganizer.ADDRBOOK,
			numTotal: null, //we don't know how many
			noTooltip: true //so don't show tooltip
		};

		dlsFolder = new ZmAddrBook(params);
		root.children.add(dlsFolder);
		dlsFolder._isDL = true;
	}
	else {
		//the dls folder object exists but no longer as a child of the root.
		dlsFolder.parent = root;
		root.children.add(dlsFolder); //any better way to do this?
	}

};

ZmContactList.prototype.add = 
function(item, index) {
	if (!item.id || !this._idHash[item.id]) {
		this._vector.add(item, index);
		if (item.id) {
			this._idHash[item.id] = item;
		}
		this._updateHashes(item, true);
	}
};

ZmContactList.prototype.cache = 
function(offset, newList) {
	var getId = function(){
		return this.id;
	}
	var exists = function(obj) {
		return this._vector.containsLike(obj, getId);
	}
	var unique = newList.sub(exists, this);

	this.getVector().merge(offset, unique);
	// reparent each item within new list, and add it to ID hash
	var list = unique.getArray();
	for (var i = 0; i < list.length; i++) {
		var item = list[i];
		item.list = this;
		if (item.id) {
			this._idHash[item.id] = item;
		}
	}
};

/**
 * @private
 */
ZmContactList.prototype._addContact =
function(contact) {

	// note that we don't create a ZmContact here (optimization)
	contact.list = this;
	this._updateHashes(contact, true);
	var fn = [], fl = [];
	if (contact._attrs[ZmContact.F_firstName])	{ fn.push(contact._attrs[ZmContact.F_firstName]); }
	if (contact._attrs[ZmContact.F_middleName])	{ fn.push(contact._attrs[ZmContact.F_middleName]); }
	if (contact._attrs[ZmContact.F_lastName])	{ fn.push(contact._attrs[ZmContact.F_lastName]); }
	if (fn.length) {
		contact._attrs[ZmContact.X_fullName] = fn.join(" ");
	}
	if (contact._attrs[ZmContact.F_firstName])	{ fl.push(contact._attrs[ZmContact.F_firstName]); }
	if (contact._attrs[ZmContact.F_lastName])	{ fl.push(contact._attrs[ZmContact.F_lastName]); }
	contact._attrs[ZmContact.X_firstLast] = fl.join(" ");

	this.add(contact);
};

/**
 * Converts an anonymous contact object (contained by the JS returned by load request)
 * into a ZmContact, and updates the containing list if it is the canonical one.
 *
 * @param {Object}	contact		a contact
 * @param {int}	idx		the index of contact in canonical list
 * 
 * @private
 */
ZmContactList.prototype._realizeContact =
function(contact, idx) {

	if (contact instanceof ZmContact) { return contact; }
	if (contact && contact.type == ZmItem.CONTACT) { return contact; }	// instanceof often fails in new window

	var args = {list:this};
	var obj = eval(ZmList.ITEM_CLASS[this.type]);
	var realContact = obj && obj.createFromDom(contact, args);

	if (this.isCanonical) {
		var a = this.getArray();
		idx = idx || this.getIndexById(contact.id);
		a[idx] = realContact;
		this._updateHashes(realContact, true);
		this._idHash[contact.id] = realContact;
	}

	return realContact;
};

/**
 * Finds the array index for the contact with the given ID.
 *
 * @param {int}	id		the contact ID
 * @return	{int}	the index
 * @private
 */
ZmContactList.prototype.getIndexById =
function(id) {
	var a = this.getArray();
	for (var i = 0; i < a.length; i++) {
		if (a[i].id == id) {
			return i;
		}
	}
	return null;
};

/**
 * Override in order to make sure the contacts have been realized. We don't
 * call realizeContact() since this is not the canonical list.
 *
 * @param {int}	offset		the starting index
 * @param {int}	limit		the size of sublist
 * @return	{AjxVector}	a vector of {@link ZmContact} objects
 */
ZmContactList.prototype.getSubList =
function(offset, limit, folderId) {
	if (folderId && this.isCanonical) {
		// only collect those contacts that belong to the given folderId if provided
		var newlist = [];
		var sublist = this.getArray();
		var offsetCount = 0;
		this.setHasMore(false);

		for (var i = 0; i < sublist.length; i++) {
			sublist[i] = this._realizeContact(sublist[i], i);
			var folder = appCtxt.getById(sublist[i].folderId);
			if (folder && folder.nId == ZmOrganizer.normalizeId(folderId)) {
				if (offsetCount >= offset) {
					if (newlist.length == limit) {
						this.setHasMore(true);
						break;
					}
					newlist.push(sublist[i]);
				}
				offsetCount++;
			}
		}

		return AjxVector.fromArray(newlist);
	} else {
		var vec = ZmList.prototype.getSubList.call(this, offset, limit);
		if (vec) {
			var a = vec.getArray();
			for (var i = 0; i < a.length; i++) {
				a[i] = this._realizeContact(a[i], offset + i);
			}
		}

		return vec;
	}
};

/**
 * Override in order to make sure the contact has been realized. Canonical list only.
 *
 * @param {int}	id		the contact ID
 * @return	{ZmContact}	the contact or <code>null</code> if not found
 */
ZmContactList.prototype.getById =
function(id) {
	if (!id || !this.isCanonical) return null;

	var contact = this._idHash[id];
	return contact ? this._realizeContact(contact) : null;
};

/**
 * Gets the contact with the given address, if any (canonical list only).
 *
 * @param {String}	address	an email address
 * @return	{ZmContact}	the contact or <code>null</code> if not found
 */
ZmContactList.prototype.getContactByEmail =
function(address) {
	if (!address || !this.isCanonical) return null;

	var contact = this._emailToContact[address.toLowerCase()];
	if (contact) {
		contact = this._realizeContact(contact);
		contact._lookupEmail = address;	// so caller knows which address matched
		return contact;
	} else {
		return null;
	}
};

/**
 * Gets information about the contact with the given phone number, if any (canonical list only).
 *
 * @param {String}	phone	the phone number
 * @return	{Hash}	an object with <code>contact</code> = the contact & <code>field</code> = the field with the matching phone number
 */
ZmContactList.prototype.getContactByPhone =
function(phone) {
	if (!phone || !this.isCanonical) return null;

	var digits = this._getPhoneDigits(phone);
	var data = this._phoneToContact[digits];
	if (data) {
		data.contact = this._realizeContact(data.contact);
		return data;
	} else {
		return null;
	}
};

/**
 * Moves a list of items to the given folder.
 * <p>
 * This method calls the base class for normal "moves" UNLESS we're dealing w/
 * shared items (or folder) in which case we must send a CREATE request for the
 * given folder to the server followed by a hard delete of the shared contact.
 * </p>
 *
 * @param {Hash}	params		a hash of parameters
 * @param	{Array}       params.items			a list of items to move
 * @param	{ZmFolder}	params.folder		the destination folder
 * @param	{Hash}	       params.attrs		the additional attrs for SOAP command
 * @param	{Boolean}	params.outOfTrash	if <code>true</code>, we are moving contacts out of trash
 */
ZmContactList.prototype.moveItems =
function(params) {

	params = Dwt.getParams(arguments, ["items", "folder", "attrs", "outOfTrash"]);
	params.items = AjxUtil.toArray(params.items);

	var moveBatchCmd = new ZmBatchCommand(true, null, true);
	var loadBatchCmd = new ZmBatchCommand(true, null, true);
	var softMove = [];

	// if the folder we're moving contacts to is a shared folder, then dont bother
	// checking whether each item is shared or not
	if (params.items[0] && params.items[0] instanceof ZmItem) {
		for (var i = 0; i < params.items.length; i++) {
			var contact = params.items[i];

			if (contact.isReadOnly()) { continue; }

			softMove.push(contact);
		}
	} else {
		softMove = params.items;
	}

	// for "soft" moves, handle moving out of Trash differently
	if (softMove.length > 0) {
		var params1 = AjxUtil.hashCopy(params);
		params1.attrs = params.attrs || {};
		var toFolder = params.folder;
		params1.attrs.l = toFolder.isRemote() ? toFolder.getRemoteId() : toFolder.id;
		params1.action = "move";
        params1.accountName = appCtxt.multiAccounts && appCtxt.accountList.mainAccount.name;
        if (params1.folder.id == ZmFolder.ID_TRASH) {
            params1.actionTextKey = 'actionTrash';
            // bug: 47389 avoid moving to local account's Trash folder.
            params1.accountName = appCtxt.multiAccounts && params.items[0].getAccount().name;
        } else {
            params1.actionTextKey = 'actionMove';
            params1.actionArg = toFolder.getName(false, false, true);
        }
		params1.callback = params.outOfTrash && new AjxCallback(this, this._handleResponseMoveItems, params);

		this._itemAction(params1);
	}
};

/**
 * @private
 */
ZmContactList.prototype._handleResponseMoveBatchCmd =
function(result) {
	var resp = result.getResponse().BatchResponse.ContactActionResponse;
	// XXX: b/c the server does not return notifications for actions done on
	//      shares, we manually notify - TEMP UNTIL WE GET BETTER SERVER SUPPORT
	var ids = resp[0].action.id.split(",");
	for (var i = 0; i < ids.length; i++) {
		var contact = appCtxt.cacheGet(ids[i]);
		if (contact && contact.isShared()) {
			contact.notifyDelete();
			appCtxt.cacheRemove(ids[i]);
		}
	}
};

/**
 * @private
 */
ZmContactList.prototype._handleResponseLoadMove =
function(moveBatchCmd, params) {
	var deleteCmd = new AjxCallback(this, this._itemAction, [params]);
	moveBatchCmd.add(deleteCmd);

	var respCallback = new AjxCallback(this, this._handleResponseMoveBatchCmd);
	moveBatchCmd.run(respCallback);
};

/**
 * @private
 */
ZmContactList.prototype._handleResponseBatchLoad =
function(batchCmd, folder, result, contact) {
	batchCmd.add(this._getCopyCmd(contact, folder));
};

/**
 * @private
 */
ZmContactList.prototype._getCopyCmd =
function(contact, folder) {
	var temp = new ZmContact(null, this);
	for (var j in contact.attr) {
		temp.attr[j] = contact.attr[j];
	}
	temp.attr[ZmContact.F_folderId] = folder.id;

	return new AjxCallback(temp, temp.create, [temp.attr]);
};

/**
 * Deletes contacts after checking that this is not a GAL list.
 *
 * @param {Hash}	params		a hash of parameters
 * @param	{Array}	       params.items			the list of items to delete
 * @param	{Boolean}	params.hardDelete	if <code>true</code>, force physical removal of items
 * @param	{Object}	params.attrs			the additional attrs for SOAP command
 */
ZmContactList.prototype.deleteItems =
function(params) {
	if (this.isGal) {
		if (ZmContactList.deleteGalItemsAllowed(params.items)) {
			this._deleteDls(params.items);
			return;
		}
		DBG.println(AjxDebug.DBG1, "Cannot delete GAL contacts that are not DLs");
		return;
	}
	ZmList.prototype.deleteItems.call(this, params);
};

ZmContactList.deleteGalItemsAllowed =
function(items) {
	var deleteDomainsAllowed = appCtxt.createDistListAllowedDomainsMap;
	if (items.length == 0) {
		return false; //need a special case since we don't want to enable the "delete" button for 0 items.
	}
	for (var i = 0; i < items.length; i++) {
		var contact = items[i];
		var email = contact.getEmail();
		var domain = email.split("@")[1];
		var isDL = contact && contact.isDistributionList();
		//see bug 71368 and also bug 79672 - the !contact.dlInfo is in case somehow dlInfo is missing - so unfortunately if that happens (can't repro) - let's not allow to delete since we do not know if it's an owner
		if (!isDL || !deleteDomainsAllowed[domain] || !contact.dlInfo || !contact.dlInfo.isOwner) {
			return false;
		}
	}
	return true;
};

ZmContactList.prototype._deleteDls =
function(items, confirmDelete) {

	if (!confirmDelete) {
		var callback = this._deleteDls.bind(this, items, true);
		this._popupDeleteWarningDialog(callback, false, items.length);
		return;
	}

	var reqs = [];
	for (var i = 0; i < items.length; i++) {
		var contact = items[i];
		var email = contact.getEmail();
		reqs.push({
				_jsns: "urn:zimbraAccount",
				dl: {by: "name",
					 _content: contact.getEmail()
				},
				action: {
					op: "delete"
				}
			});
	}
	var jsonObj = {
		BatchRequest: {
			_jsns: "urn:zimbra",
			DistributionListActionRequest: reqs
		}
	};
	var respCallback = this._deleteDlsResponseHandler.bind(this, items);
	appCtxt.getAppController().sendRequest({jsonObj: jsonObj, asyncMode: true, callback: respCallback});

};

ZmContactList.prototype._deleteDlsResponseHandler =
function(items) {
	if (appCtxt.getCurrentView().isZmGroupView) {
		//this is the case we were editing the DL (different than viewing it in the DL list, in which case it's the contactListController).
		//so we now need to pop up the view.
		this.controller.popView();
	}

	appCtxt.setStatusMsg(items.length == 1 ? ZmMsg.dlDeleted : ZmMsg.dlsDeleted);

	for (var i = 0; i < items.length; i++) {
		var item = items[i];
		item.clearDlInfo();
		item._notify(ZmEvent.E_DELETE);
	}
};



/**
 * Sets the is GAL flag.
 * 
 * @param	{Boolean}	isGal		<code>true</code> if contact list is GAL
 */
ZmContactList.prototype.setIsGal =
function(isGal) {
	this.isGal = isGal;
};

ZmContactList.prototype.notifyCreate =
function(node) {
	var obj = eval(ZmList.ITEM_CLASS[this.type]);
	if (obj) {
		var item = obj.createFromDom(node, {list:this});
		var index = this._sortIndex(item);
		// only add if it sorts into this list
		var listSize = this.size();
		var visible = false;
		if (index < listSize || listSize == 0 || (index==listSize && !this._hasMore)) {
			this.add(item, index);
			this.createLocal(item);
			visible = true;
		}
		this._notify(ZmEvent.E_CREATE, {items: [item], sortIndex: index, visible: visible});
	}
};

/**
 * Moves the items.
 * 
 * @param	{Array}	items		an array of {@link ZmContact} objects
 * @param	{String}	folderId	the folder id
 */
ZmContactList.prototype.moveLocal =
function(items, folderId) {
	// don't remove any contacts from the canonical list
	if (!this.isCanonical)
		ZmList.prototype.moveLocal.call(this, items, folderId);
	if (folderId == ZmFolder.ID_TRASH) {
		for (var i = 0; i < items.length; i++) {
			this._updateHashes(items[i], false);
		}
	}
};

/**
 * Deletes the items.
 * 
 * @param	{Array}	items		an array of {@link ZmContact} objects
 */
ZmContactList.prototype.deleteLocal =
function(items) {
	ZmList.prototype.deleteLocal.call(this, items);
	for (var i = 0; i < items.length; i++) {
		this._updateHashes(items[i], false);
	}
};

/**
 * Handle modified contact.
 * 
 * @private
 */
ZmContactList.prototype.modifyLocal =
function(item, details) {
	if (details) {
		// notify item's list
		this._evt.items = details.items = [item];
		this._evt.item = details.contact; //somehow this was set to something obsolete. What a mess. Also note that item is Object while details.contact is ZmContact
		this._notify(ZmEvent.E_MODIFY, details);
	}

	var contact = details.contact;
	if (this.isCanonical || contact.attr[ZmContact.F_email] != details.oldAttr[ZmContact.F_email]) {
		// Remove traces of old contact - NOTE: we pass in null for the ID on
		// PURPOSE to avoid overwriting the existing cached contact
		var oldContact = new ZmContact(null, this);
		oldContact.id = details.contact.id;
		oldContact.attr = details.oldAttr;
		this._updateHashes(oldContact, false);

		// add new contact to hashes
		this._updateHashes(contact, true);
	}

	// place in correct position in list
	if (details.fileAsChanged) {
		this.remove(contact);
		var index = this._sortIndex(contact);
		var listSize = this.size();
		if (index < listSize || listSize == 0 || (index == listSize && !this._hasMore)) {
			this.add(contact, index);
		}
	}

	// reset addrbook property
	if (contact.addrbook && (contact.addrbook.id != contact.folderId)) {
		contact.addrbook = appCtxt.getById(contact.folderId);
	}
};

/**
 * Creates the item local.
 * 
 * @param	{ZmContact}	item		the item
 */
ZmContactList.prototype.createLocal =
function(item) {
	this._updateHashes(item, true);
};

/**
 * @private
 */
ZmContactList.prototype._updateHashes =
function(contact, doAdd) {

	this._app.updateCache(contact, doAdd);

	// Update email hash.
	for (var index = 0; index < ZmContact.EMAIL_FIELDS.length; index++) {
		var field = ZmContact.EMAIL_FIELDS[index];
		for (var i = 1; true; i++) {
			var aname = ZmContact.getAttributeName(field, i);
			var avalue = ZmContact.getAttr(contact, aname);
			if (!avalue) break;
			if (doAdd) {
				this._emailToContact[avalue.toLowerCase()] = contact;
			} else {
				delete this._emailToContact[avalue.toLowerCase()];
			}
		}
	}

	// Update phone hash.
	if (appCtxt.get(ZmSetting.VOICE_ENABLED) || this._alwaysUpdateHashes) {
		for (var index = 0; index < ZmContact.PHONE_FIELDS.length; index++) {
			var field = ZmContact.PHONE_FIELDS[index];
			for (var i = 1; true; i++) {
				var aname = ZmContact.getAttributeName(field, i);
				var avalue = ZmContact.getAttr(contact, aname);
				if (!avalue) break;
				var digits = this._getPhoneDigits(avalue);
				if (digits) {
					if (doAdd) {
						this._phoneToContact[avalue] = {contact: contact, field: aname};
					} else {
						delete this._phoneToContact[avalue];
					}
				}
			}
		}
	}
};

/**
 * Strips all non-digit characters from a phone number.
 * 
 * @private
 */
ZmContactList.prototype._getPhoneDigits =
function(phone) {
	return phone.replace(/[^\d]/g, '');
};

/**
 * Returns the position at which the given contact should be inserted in this list.
 * 
 * @private
 */
ZmContactList.prototype._sortIndex =
function(contact) {
	var a = this._vector.getArray();
	for (var i = 0; i < a.length; i++) {
		if (ZmContact.compareByFileAs(a[i], contact) > 0) {
			return i;
		}
	}
	return a.length;
};

/**
 * Gets the list ID hash
 * @return idHash {Ojbect} list ID hash
 */
ZmContactList.prototype.getIdHash =
function() {
	return this._idHash;
}

/**
 * @private
 */
ZmContactList.prototype._handleResponseModifyItem =
function(item, result) {
	// NOTE: we overload and do nothing b/c base class does more than we want
	//       (since everything is handled by notifications)
};
}
if (AjxPackage.define("zimbraMail.abook.view.ZmContactsHelper")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */

/**
 * @overview
 * This file contains a contact helper class.
 * 
 */

/**
 * Default constructor for helper class.
 * @class
 * Miscellaneous contacts-related utility functions. So far, mostly things that
 * {@link ZmContactPicker} and {@link ZmGroupView} both need to perform a contacts search and
 * display results in a list view.
 *
 * @author Conrad Damon
 */
ZmContactsHelper = function() {};

/**
 * Performs a contact search (in either personal contacts or in the GAL) and populates
 * the source list view with the results.
 *
 * @param	{Hash}	params		a hash of parameters
 * @param {Object}	params.obj			the object that is doing the search
 * @param {String}	params.query			the query string to search on
 * @param {String}	params.queryHint		the query hint (i.e. searching shared folders)
 * @param {Boolean}	params.ascending		if <code>true</code>, sort in ascending order
 * @param {int}	params.lastId		the ID of last item displayed (for pagination)
 * @param {String}	params.lastSortVal	the value of sort field for above item
 * @param {AjxCallback}	params.respCallback	the callback to call once response comes back from server
 * @param {AjxCallback}	params.errorCallback	the callback to call if error returned from server
 * @param {String}	params.accountName	the account to make search request on behalf of
 * @param {Array}   params.conds	the conds to restrict the search by (array of {attr:"", op:"", value:""} hashes)
 */
ZmContactsHelper.search =
function(params) {
	var o = params.obj;
	if (o._searchButton) {
		o._searchButton.setEnabled(false);
	}

	params.sortBy = params.ascending ? ZmSearch.NAME_ASC : ZmSearch.NAME_DESC;
	params.types = AjxVector.fromArray([ZmItem.CONTACT]);
	params.offset = params.offset || 0;
	params.limit = ZmContactsApp.SEARCHFOR_MAX;
	params.contactSource = o._contactSource;
	params.field = "contact";

	var search = new ZmSearch(params);
	search.execute({callback:params.respCallback, errorCallback:params.errorCallback});
};

/**
 * Take the contacts and create a list of their email addresses (a contact may have more than one)
 * 
 * @private
 */
ZmContactsHelper._processSearchResponse = 
function(resp, includeContactsWithNoEmail) {
	var vec = resp.getResults(ZmItem.CONTACT);

	// Take the contacts and create a list of their email addresses (a contact may have more than one)
	var list = [];
	var a = vec.getArray();
	for (var i = 0; i < a.length; i++) {
		var contact = a[i];
		if (contact.isGroup() && !contact.isDL) {
			var members = contact.getGroupMembers().good.toString(AjxEmailAddress.SEPARATOR);
			ZmContactsHelper._addContactToList(list, contact, members, true);
		} else {
			var emails = contact.isGal ? [contact.getEmail()] : contact.getEmails();
			for (var j = 0; j < emails.length; j++) {
				ZmContactsHelper._addContactToList(list, contact, emails[j]);
			}
			if (includeContactsWithNoEmail && emails.length == 0) {
				ZmContactsHelper._addContactToList(list, contact, null);
			}
		}
	}
	
	return list;
};

/**
 * @private
 */
ZmContactsHelper._addContactToList = 
function(list, contact, addr, isGroup) {

	var email = ZmContactsHelper._wrapContact(contact, addr, isGroup);  
	list.push(email);
};

/**
 * wrapps the contact inside a AjxEmailAddress object, and adds a couple extra fields to the AjxEmailAddress instance (value, contact, icon [which I'm not sure is used])
 *
 * @param contact
 * @param addr {String} optional.
 * @param isGroup
 */
ZmContactsHelper._wrapContact =
function(contact, addr, isGroup) {

	addr = addr || contact.getEmail();
	var fileAs = contact.getFileAs();
	var name = (fileAs != addr) ? fileAs : "";  //todo ??? this is weird.
	var	type = contact.isGal ? ZmContact.GROUP_GAL_REF : ZmContact.GROUP_CONTACT_REF;
	var	value = contact.isGal ? (contact.ref || contact.id) : contact.id;  //defaulting to contact.id in the gal case since from GetContactsResponse the ref is not returned and we can end up with it cached without the ref. Probably need to fix that.
	var displayName = contact.getFullNameForDisplay();

	var email = new AjxEmailAddress(addr, type, name, displayName, isGroup);

	email.value = value;
	email.id = Dwt.getNextId();
	email.__contact = contact;
	email.icon = contact.getIcon();
	if (contact.isDL) {
		email.isGroup = true;
		email.canExpand = contact.canExpand;
		var ac = window.parentAppCtxt || window.appCtxt;
		ac.setIsExpandableDL(addr, email.canExpand);
	}
	return email;
};

/**
 * wrapps the inline address (there's no real ZmContact object) inside AjxEmailAddress and adds the value attribute to it.
 * this is so we treat real contacts and inline contacts consistently throughout the rest of the code.
 *
 * @param value  {String} - the inline email address and/or name (e.g. "john doe <john@doe.com>" or "john@doe.com")
 */
ZmContactsHelper._wrapInlineContact =
function(value) {
	var email = AjxEmailAddress.parse(value); //from legacy data at least (not sure about new), the format might be something like "Inigo Montoya <inigo@theprincessbride.com>" so we have to parse.
	if (!email) {
		//this can happen when creating inline in contact group edit, and the user did not suply email address in the inline value
		email = new AjxEmailAddress(value, null, value);
	}
	email.type = ZmContact.GROUP_INLINE_REF;
	email.value = value;
	email.id = Dwt.getNextId();
	return email;
};


/**
 * The items are AjxEmailAddress objects
 * 
 * @private
 */
ZmContactsHelper._getEmailField =
function(html, idx, item, field, colIdx) {
	if (field == ZmItem.F_TYPE) {
		html[idx++] = AjxImg.getImageHtml(item.icon);
	} else if (field == ZmItem.F_NAME) {
		html[idx++] = '<span style="white-space:nowrap">';
		html[idx++] = AjxStringUtil.htmlEncode(item.name || ZmMsg.noName);
		html[idx++] = "</span>";
	} else if (field == ZmItem.F_EMAIL) {
		html[idx++] = AjxStringUtil.htmlEncode(item.address);
	} else if (field == ZmItem.F_DEPARTMENT) {
		if (item.__contact) {
			html[idx++] = AjxStringUtil.htmlEncode(ZmContact.getAttr(item.__contact, ZmContact.F_department));
		}
	}
	return idx;
};
}
if (AjxPackage.define("zimbraMail.abook.view.ZmContactPicker")) {
/*
 * ***** BEGIN LICENSE BLOCK *****
 * Zimbra Collaboration Suite Web Client
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc.
 *
 * The contents of this file are subject to the Common Public Attribution License Version 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at: https://www.zimbra.com/license
 * The License is based on the Mozilla Public License Version 1.1 but Sections 14 and 15
 * have been added to cover use of software over a computer network and provide for limited attribution
 * for the Original Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied.
 * See the License for the specific language governing rights and limitations under the License.
 * The Original Code is Zimbra Open Source Web Client.
 * The Initial Developer of the Original Code is Zimbra, Inc.  All rights to the Original Code were
 * transferred by Zimbra, Inc. to Synacor, Inc. on September 14, 2015.
 *
 * All portions of the code are Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016 Synacor, Inc. All Rights Reserved.
 * ***** END LICENSE BLOCK *****
 */

/**
 * @overview
 * This file contains the contact picker classes.
 * 
 */

/**
 * Creates a dialog that lets the user select addresses from a contact list.
 * @constructor
 * @class
 * This class creates and manages a dialog that lets the user select addresses
 * from a contact list. Two lists are maintained, one with contacts to select
 * from, and one that contains the selected addresses. Between them are buttons
 * to shuffle addresses back and forth between the two lists.
 *
 * @author Conrad Damon
 * 
 * @param {Array}	buttonInfo		the transfer button IDs and labels
 * 
 * @extends		DwtDialog
 */
ZmContactPicker = function(buttonInfo) {

	DwtDialog.call(this, {parent:appCtxt.getShell(), title:ZmMsg.selectAddresses, id: "ZmContactPicker"});

	this._buttonInfo = buttonInfo;
	this._initialized = false;
	this._emailListOffset = 0; //client side paginating over email list. Offset of current page of email addresses. Quite different than _lastServerOffset if contacts have 0 or more than 1 email addresses.
	this._serverContactOffset = 0; //server side paginating over contact list. Offset of last contact block we got from the server (each contact could have 0, 1, or more emails so we have to keep track of this separate from the view offset.
	this._ascending = true; //asending or descending search. Keep it stored for pagination to do the right sort.
	this._emailList = new AjxVector();
	this._detailedSearch = appCtxt.get(ZmSetting.DETAILED_CONTACT_SEARCH_ENABLED);
	this._ignoreSetDragBoundries = true;

	this.setSize(Dwt.DEFAULT, this._getDialogHeight());

	this._searchErrorCallback = new AjxCallback(this, this._handleErrorSearch);
};

ZmContactPicker.prototype = new DwtDialog;
ZmContactPicker.prototype.constructor = ZmContactPicker;

ZmContactPicker.prototype.isZmContactPicker = true;
ZmContactPicker.prototype.toString = function() { return "ZmContactPicker"; };

// Consts

ZmContactPicker.DIALOG_HEIGHT = 460;

ZmContactPicker.SEARCH_BASIC = "search";
ZmContactPicker.SEARCH_NAME = "name";
ZmContactPicker.SEARCH_EMAIL = "email";
ZmContactPicker.SEARCH_DEPT = "dept";
ZmContactPicker.SEARCH_PHONETIC = "phonetic";

ZmContactPicker.SHOW_ON_GAL = [ZmContactPicker.SEARCH_BASIC, ZmContactPicker.SEARCH_NAME, ZmContactPicker.SEARCH_EMAIL, ZmContactPicker.SEARCH_DEPT];
ZmContactPicker.SHOW_ON_NONGAL = [ZmContactPicker.SEARCH_BASIC, ZmContactPicker.SEARCH_NAME, ZmContactPicker.SEARCH_PHONETIC, ZmContactPicker.SEARCH_EMAIL];
ZmContactPicker.ALL = [ ZmContactPicker.SEARCH_BASIC, ZmContactPicker.SEARCH_NAME, ZmContactPicker.SEARCH_PHONETIC, ZmContactPicker.SEARCH_EMAIL, ZmContactPicker.SEARCH_DEPT ];

// Public methods


/**
* Displays the contact picker dialog. The source list is populated with
* contacts, and the target list is populated with any addresses that are
* passed in. The address button that was used to popup the dialog is set
* as the active button.
*
* @param {String}	buttonId	the button ID of the button that called us
* @param {Hash}	addrs		a hash of 3 vectors (one for each type of address)
* @param {String}	str		initial search string
*/
ZmContactPicker.prototype.popup =
function(buttonId, addrs, str, account) {

	if (!this._initialized) {
		this._initialize(account);
		this._initialized = true;
	}
	else if (appCtxt.multiAccounts && this._account != account) {
		this._account = account;
		this._resetSelectDiv();
	}
	this._emailListOffset = 0;

	var searchFor = this._searchInSelect ? this._searchInSelect.getValue() : ZmContactsApp.SEARCHFOR_CONTACTS;

	// reset column sorting preference
	this._chooser.sourceListView.setSortByAsc(ZmItem.F_NAME, true);

	// reset button states
	this._chooser.reset();
	if (buttonId) {
		this._chooser._setActiveButton(buttonId);
	}

	// populate target list if addrs were passed in
	if (addrs) {
		for (var id in addrs) {
			this._chooser.addItems(addrs[id], DwtChooserListView.TARGET, true, id);
		}
	}

	for (var fieldId in this._searchField) {
		var field = this._searchField[fieldId];
		field.disabled = false;
		field.value = (AjxUtil.isObject(str) ? str[fieldId] : str) || "";
	}

	// reset paging buttons
	this._prevButton.setEnabled(false);
	this._nextButton.setEnabled(false);

	this.search(null, true, true);

    DwtDialog.prototype.popup.call(this);

	this._resizeChooser();

    if ((this.getLocation().x < 0 ||  this.getLocation().y < 0) ){
                // parent window size is smaller than Dialog size
                this.setLocation(0,30);
                var size = Dwt.getWindowSize();
                var currentSize = this.getSize();
                var dragElement = document.getElementById(this._dragHandleId);
                DwtDraggable.setDragBoundaries(dragElement, 100 - currentSize.x, size.x - 100, 0, size.y - 100);
    }

	var focusField = this._searchField[ZmContactPicker.SEARCH_BASIC] || this._searchField[ZmContactPicker.SEARCH_NAME];
	appCtxt.getKeyboardMgr().grabFocus(focusField);

};


ZmContactPicker.prototype._resetResults =
function() {
	this._emailList.removeAll();
	this._serverContactOffset = 0;
	this._emailListOffset = 0;
};

/**
 * Closes the dialog.
 * 
 */
ZmContactPicker.prototype.popdown =
function() {
	// disable search field (hack to fix bleeding cursor)

	for (var fieldId in this._searchField) {
		this._searchField[fieldId].disabled = true;
	}

	this._contactSource = null;
	this._resetResults();

	DwtDialog.prototype.popdown.call(this);
};

/**
 * Performs a search.
 * 
 * @private
 */
ZmContactPicker.prototype.search =
function(colItem, ascending, firstTime, lastId, lastSortVal, offset) {
	if (offset == undefined) {
		//this could be a call from DwtChooserListView.prototype._sortColumn, which means we have to reset the result and both server and client pagination.
		//In any case the results should be reset or are already reset so doesn't hurt to reset.
		this._resetResults();
	}

	if (ascending === null || ascending === undefined) {
		ascending = this._ascending;
	}
	else {
		this._ascending = ascending;
	}
	
	var query;
	var queryHint = [];
	var emailQueryTerm = "";
	var phoneticQueryTerms = [];
	var nameQueryTerms = [];
	var conds = [];
	if (this._detailedSearch) {
		var nameQuery = this.getSearchFieldValue(ZmContactPicker.SEARCH_NAME);
		var emailQuery = this.getSearchFieldValue(ZmContactPicker.SEARCH_EMAIL);
		var deptQuery = this.getSearchFieldValue(ZmContactPicker.SEARCH_DEPT);
		var phoneticQuery = this.getSearchFieldValue(ZmContactPicker.SEARCH_PHONETIC);
		var isGal = this._searchInSelect && (this._searchInSelect.getValue() == ZmContactsApp.SEARCHFOR_GAL);
		if (nameQuery && !isGal) {
			var nameQueryPieces = nameQuery.split(/\s+/);
			for (var i = 0; i < nameQueryPieces.length; i++) {
				var nameQueryPiece = nameQueryPieces[i];
				nameQueryTerms.push("#"+ZmContact.F_firstName + ":" + nameQueryPiece);
				nameQueryTerms.push("#"+ZmContact.F_lastName + ":" + nameQueryPiece);
				nameQueryTerms.push("#"+ZmContact.F_middleName + ":" + nameQueryPiece);
				nameQueryTerms.push("#"+ZmContact.F_nickname + ":" + nameQueryPiece);
			}
			query = "(" + nameQueryTerms.join(" OR ") + ")";
		} else {
			if (nameQuery && isGal) {
				conds.push([{attr:ZmContact.F_firstName, op:"has", value: nameQuery},
				            {attr:ZmContact.F_lastName,  op:"has", value: nameQuery},
				            {attr:ZmContact.F_middleName, op:"has", value: nameQuery},
				            {attr:ZmContact.F_nickname,  op:"has", value: nameQuery},
				            {attr:ZmContact.F_phoneticFirstName, op:"has", value: nameQuery},
				            {attr:ZmContact.F_phoneticLastName,  op:"has", value: nameQuery}]);
			}
			query = "";
		}
		if (emailQuery) {
			if (isGal) {
				conds.push([{attr:ZmContact.F_email, op:"has", value: emailQuery},
				{attr:ZmContact.F_email2, op:"has", value: emailQuery},
				{attr:ZmContact.F_email3, op:"has", value: emailQuery},
				{attr:ZmContact.F_email4, op:"has", value: emailQuery},
				{attr:ZmContact.F_email5, op:"has", value: emailQuery},
				{attr:ZmContact.F_email6, op:"has", value: emailQuery},
				{attr:ZmContact.F_email7, op:"has", value: emailQuery},
				{attr:ZmContact.F_email8, op:"has", value: emailQuery},
				{attr:ZmContact.F_email9, op:"has", value: emailQuery},
				{attr:ZmContact.F_email10, op:"has", value: emailQuery},
				{attr:ZmContact.F_email11, op:"has", value: emailQuery},
				{attr:ZmContact.F_email12, op:"has", value: emailQuery},
				{attr:ZmContact.F_email13, op:"has", value: emailQuery},
				{attr:ZmContact.F_email14, op:"has", value: emailQuery},
				{attr:ZmContact.F_email15, op:"has", value: emailQuery},
				{attr:ZmContact.F_email16, op:"has", value: emailQuery}
				]);
			} else {
				emailQueryTerm = "to:"+emailQuery+"*";
			}
		}
		if (deptQuery && isGal) {
			conds.push({attr:ZmContact.F_department, op:"has", value: deptQuery});
		}
		if (phoneticQuery && !isGal) {
			var phoneticQueryPieces = phoneticQuery.split(/\s+/);
			for (var i=0; i<phoneticQueryPieces.length; i++) {
				phoneticQueryTerms.push("#"+ZmContact.F_phoneticFirstName + ":" + phoneticQueryPieces[i]);
				phoneticQueryTerms.push("#"+ZmContact.F_phoneticLastName + ":" + phoneticQueryPieces[i]);
			}
		}
	} else {
		query = this.getSearchFieldValue(ZmContactPicker.SEARCH_BASIC);
	}
	

	if (this._searchInSelect) {
		var searchFor = this._searchInSelect.getValue();
		this._contactSource = (searchFor == ZmContactsApp.SEARCHFOR_CONTACTS || searchFor == ZmContactsApp.SEARCHFOR_PAS)
			? ZmItem.CONTACT
			: ZmId.SEARCH_GAL;

		if (searchFor == ZmContactsApp.SEARCHFOR_PAS) {
			queryHint.push(ZmSearchController.generateQueryForShares(ZmId.ITEM_CONTACT) || "is:local");
		} else if (searchFor == ZmContactsApp.SEARCHFOR_CONTACTS) {
			queryHint.push("is:local");
		} else if (searchFor == ZmContactsApp.SEARCHFOR_GAL) {
            ascending = true;
        }
	} else {
		this._contactSource = appCtxt.get(ZmSetting.CONTACTS_ENABLED, null, this._account)
			? ZmItem.CONTACT
			: ZmId.SEARCH_GAL;

		if (this._contactSource == ZmItem.CONTACT) {
			queryHint.push("is:local");
		}
	}

	if (this._contactSource == ZmItem.CONTACT && query != "" && !query.startsWith ("(")) {
		query = query.replace(/\"/g, '\\"');
		query = "\"" + query + "\"";
	}
	if (phoneticQueryTerms.length) {
		query = query + " (" + phoneticQueryTerms.join(" OR ") + ")";
	}
	if (emailQueryTerm.length) {
		query = query + " " + emailQueryTerm;  // MUST match email term, hence AND rather than OR
	}

	if (this._searchIcon) { //does not exist in ZmGroupView case
		this._searchIcon.className = "DwtWait16Icon";
	}

	// XXX: line below doesn't have intended effect (turn off column sorting for GAL search)
	if (this._chooser) { //_chooser not defined in ZmGroupView but we also do not support sorting there anyway
		this._chooser.sourceListView.sortingEnabled = (this._contactSource == ZmItem.CONTACT);
	}

	var params = {
		obj:			this,
		ascending:		ascending,
		query:			query,
		queryHint:		queryHint.join(" "),
		conds:			conds,
		offset:			offset || 0,
		lastId:			lastId,
		lastSortVal:	lastSortVal,
		respCallback:	(new AjxCallback(this, this._handleResponseSearch, [firstTime])),
		errorCallback:	this._searchErrorCallback,
		accountName:	(this._account && this._account.name),
		expandDL:		true
	};
	ZmContactsHelper.search(params);
};

/**
 * @private
 */
ZmContactPicker.prototype._contentHtml =
function(account) {
	var showSelect;
	if (appCtxt.multiAccounts) {
		var list = appCtxt.accountList.visibleAccounts;
		for (var i = 0; i < list.length; i++) {
			var account = list[i];
			if (appCtxt.get(ZmSetting.CONTACTS_ENABLED, null, account) &&
				(appCtxt.get(ZmSetting.GAL_ENABLED, null, account) ||
				 appCtxt.get(ZmSetting.SHARING_ENABLED, null, account)))
			{
				showSelect = true;
				break;
			}
		}
	} else {
		showSelect = (appCtxt.get(ZmSetting.CONTACTS_ENABLED) &&
					  (appCtxt.get(ZmSetting.GAL_ENABLED) ||
					   appCtxt.get(ZmSetting.SHARING_ENABLED)));
	}

	var subs = {
		id: this._htmlElId,
		showSelect: showSelect,
		detailed: this._detailedSearch
	};

	return (AjxTemplate.expand("abook.Contacts#ZmContactPicker", subs));
};

/**
 * @private
 */
ZmContactPicker.prototype._resetSelectDiv =
function() {
    this._searchInSelect.clearOptions();

    if (appCtxt.multiAccounts) {
        var accts = appCtxt.accountList.visibleAccounts;
        var org = ZmOrganizer.ITEM_ORGANIZER;
        org = ZmOrganizer.ITEM_ORGANIZER[ZmItem.CONTACT];

        for (var i = 0; i < accts.length; i++) {
            this._searchInSelect.addOption(accts[i].displayName, false, accts[i].id);
            var folderTree = appCtxt.getFolderTree(accts[i]);
            var data = [];
            data = data.concat(folderTree.getByType(org));
            for (var j = 0; j < data.length; j++) {
                var addrsbk = data[j];
                if(addrsbk.noSuchFolder) { continue; }
                this._searchInSelect.addOption(addrsbk.getName(), false, addrsbk.id, "ImgContact");
            }
            if(accts[i].isZimbraAccount && !accts[i].isMain) {
                if (appCtxt.get(ZmSetting.CONTACTS_ENABLED, null, this._account)) {
                    if (appCtxt.get(ZmSetting.SHARING_ENABLED, null, this._account))
                        this._searchInSelect.addOption(ZmMsg.searchPersonalSharedContacts, false, ZmContactsApp.SEARCHFOR_PAS, "ImgContact");
                }

                if (appCtxt.get(ZmSetting.GAL_ENABLED, null, this._account)) {
                    this._searchInSelect.addOption(ZmMsg.GAL, true, ZmContactsApp.SEARCHFOR_GAL, "ImgContact");
                }

                if (!appCtxt.get(ZmSetting.INITIALLY_SEARCH_GAL, null, this._account) ||
                        !appCtxt.get(ZmSetting.GAL_ENABLED, null, this._account))
                {
                    this._searchInSelect.setSelectedValue(ZmContactsApp.SEARCHFOR_CONTACTS);
                }
            }
        }

        for (var k = 0; k < accts.length; k++) {
            this._searchInSelect.enableOption(accts[k].id, false);
        }
    } else {

        if (appCtxt.get(ZmSetting.CONTACTS_ENABLED, null, this._account)) {
            this._searchInSelect.addOption(ZmMsg.contacts, false, ZmContactsApp.SEARCHFOR_CONTACTS);

            if (appCtxt.get(ZmSetting.SHARING_ENABLED, null, this._account))
                this._searchInSelect.addOption(ZmMsg.searchPersonalSharedContacts, false, ZmContactsApp.SEARCHFOR_PAS);
        }

        if (appCtxt.get(ZmSetting.GAL_ENABLED, null, this._account)) {
            this._searchInSelect.addOption(ZmMsg.GAL, true, ZmContactsApp.SEARCHFOR_GAL);
        }

        if (!appCtxt.get(ZmSetting.INITIALLY_SEARCH_GAL, null, this._account) ||
                !appCtxt.get(ZmSetting.GAL_ENABLED, null, this._account))
        {
            this._searchInSelect.setSelectedValue(ZmContactsApp.SEARCHFOR_CONTACTS);
        }

    }
};

ZmContactPicker.prototype.getSearchFieldValue =
function(fieldId) {
	if (!fieldId && !this._detailedSearch) {
		fieldId = ZmContactPicker.SEARCH_BASIC;
	}
	var field = this._searchField[fieldId];
	return field && AjxStringUtil.trim(field.value) || "";
};


ZmContactPicker.prototype._getDialogHeight =
function() {
	return ZmContactPicker.DIALOG_HEIGHT - (appCtxt.isChildWindow ? 100 : 0);
};

ZmContactPicker.prototype._getSectionHeight =
function(idSuffix) {
	return Dwt.getSize(document.getElementById(this._htmlElId + idSuffix)).y;

};

ZmContactPicker.prototype._resizeChooser =
function() {

	var chooserHeight = this._getDialogHeight()
			- this._getSectionHeight("_handle")  //the header
			- this._getSectionHeight("_searchTable")
			- this._getSectionHeight("_paging")
			- this._getSectionHeight("_buttonsSep")
			- this._getSectionHeight("_buttons")
			- 30; //still need some magic to account for some margins etc.

	this._chooser.resize(this.getSize().x - 25, chooserHeight);
};

/**
 * called only when ZmContactPicker is first created. Sets up initial layout.
 * 
 * @private
 */
ZmContactPicker.prototype._initialize =
function(account) {

	// create static content and append to dialog parent
	this.setContent(this._contentHtml(account));

	this._searchIcon = document.getElementById(this._htmlElId + "_searchIcon");

	// add search button
	this._searchButton = new DwtButton({parent:this, parentElement:(this._htmlElId+"_searchButton")});
	this._searchButton.setText(ZmMsg.search);
	this._searchButton.addSelectionListener(new AjxListener(this, this._searchButtonListener));

	// add select menu
	var selectCellId = this._htmlElId + "_listSelect";
	var selectCell = document.getElementById(selectCellId);
	if (selectCell) {
		this._searchInSelect = new DwtSelect({
			parent:         this,
			parentElement:  selectCellId,
			id:             Dwt.getNextId("ZmContactPickerSelect_"),
			legendId:       this._htmlElId + '_listSelectLbl'
		});
		this._resetSelectDiv();
		this._searchInSelect.addChangeListener(new AjxListener(this, this._searchTypeListener));
	} else {
		this.setSize("600");
	}

	// add chooser
	this._chooser = new ZmContactChooser({parent:this, buttonInfo:this._buttonInfo});
	this._chooser.reparentHtmlElement(this._htmlElId + "_chooser");

	// add paging buttons
	var pageListener = new AjxListener(this, this._pageListener);
	this._prevButton = new DwtButton({parent:this, parentElement:(this._htmlElId+"_pageLeft")});
	this._prevButton.setText(ZmMsg.previous);
	this._prevButton.setImage("LeftArrow");
	this._prevButton.addSelectionListener(pageListener);

	this._nextButton = new DwtButton({parent:this, style:DwtLabel.IMAGE_RIGHT, parentElement:(this._htmlElId+"_pageRight")});
	this._nextButton.setText(ZmMsg.next);
	this._nextButton.setImage("RightArrow");
	this._nextButton.addSelectionListener(pageListener);

	var pageContainer = document.getElementById(this._htmlElId + "_paging");
	if (pageContainer) {
		Dwt.setSize(pageContainer, this._chooser.sourceListView.getSize().x);
	}

	// init listeners
	this.setButtonListener(DwtDialog.OK_BUTTON, new AjxListener(this, this._okButtonListener));
	this.setButtonListener(DwtDialog.CANCEL_BUTTON, new AjxListener(this, this._cancelButtonListener));

	var fieldMap = {};
	var rowMap = {};
	this.mapFields(fieldMap, rowMap);

	this._searchField = {};
	for (var fieldId in fieldMap) {
		var field = Dwt.byId(fieldMap[fieldId]);
		if (field) {
			this._searchField[fieldId] = field;
			Dwt.setHandler(field, DwtEvent.ONKEYUP, ZmContactPicker._keyPressHdlr);
		}
	}

	this._searchRow = {};
	for (var rowId in rowMap) {
		var row = Dwt.byId(rowMap[rowId]);
		if (row) {
			this._searchRow[rowId] = row;
		}
	}
	this._updateSearchRows(this._searchInSelect && this._searchInSelect.getValue() || ZmContactsApp.SEARCHFOR_CONTACTS);
	this._keyPressCallback = new AjxCallback(this, this._searchButtonListener);
    this.sharedContactGroups = [];

	//add tabgroups for keyboard navigation
	this._tabGroup = new DwtTabGroup(this.toString());
	this._tabGroup.removeAllMembers();
	for (var i = 0; i < ZmContactPicker.ALL.length; i++) {
		field = Dwt.byId(fieldMap[ZmContactPicker.ALL[i]]);
		if (Dwt.getVisible(field)) {
			this._tabGroup.addMember(field);
		}
	}
	this._tabGroup.addMember(this._searchButton);
	this._tabGroup.addMember(this._searchInSelect);
	this._tabGroup.addMember(this._chooser.getTabGroupMember());
	this._tabGroup.addMember(this._prevButton);
	this._tabGroup.addMember(this._nextButton);
	for (var i = 0; i < this._buttonList.length; i++) {
		this._tabGroup.addMember(this._button[this._buttonList[i]]);
	}
};

ZmContactPicker.prototype.mapFields =
function(fieldMap, rowMap) {
	if (this._detailedSearch) {
		fieldMap[ZmContactPicker.SEARCH_NAME] = this._htmlElId + "_searchNameField";
		fieldMap[ZmContactPicker.SEARCH_EMAIL] = this._htmlElId + "_searchEmailField";
		fieldMap[ZmContactPicker.SEARCH_DEPT] = this._htmlElId + "_searchDepartmentField";
		fieldMap[ZmContactPicker.SEARCH_PHONETIC] = this._htmlElId + "_searchPhoneticField";
		rowMap[ZmContactPicker.SEARCH_NAME] = this._htmlElId + "_searchNameRow";
		rowMap[ZmContactPicker.SEARCH_PHONETIC] = this._htmlElId + "_searchPhoneticRow";
		rowMap[ZmContactPicker.SEARCH_EMAIL] = this._htmlElId + "_searchEmailRow";
		rowMap[ZmContactPicker.SEARCH_DEPT] = this._htmlElId + "_searchDepartmentRow";
	}
	else {
		fieldMap[ZmContactPicker.SEARCH_BASIC] = this._htmlElId + "_searchField";
		rowMap[ZmContactPicker.SEARCH_BASIC] = this._htmlElId + "_searchRow";
	}
};		

// Listeners

/**
 * @private
 */
ZmContactPicker.prototype._searchButtonListener =
function(ev) {
	this._resetResults();
	this.search();
};

/**
 * @private
 */
ZmContactPicker.prototype._handleResponseSearch =
function(firstTime, result) {
	var resp = result.getResponse();
	var serverHasMore = resp.getAttribute("more");
	var serverPaginationSupported = resp.getAttribute("paginationSupported") !== false; //if it's not specified (such as the case of SearchResponse, i.e. not Gal) it IS supported.
	this._serverHasMoreAndPaginationSupported = serverHasMore && serverPaginationSupported;
	var offset = resp.getAttribute("offset");
	this._serverContactOffset = offset || 0;
	var info = resp.getAttribute("info");
	var expanded = info && info[0].wildcard[0].expanded == "0";

	//the check for firstTime is so when the picker is popped up we probably don't want to overwhelm them with a warning message. So only show it if the user plays with the picker, using the drop-down or the search box.
	if (!firstTime && !serverPaginationSupported && (serverHasMore || expanded)) { //no idea what the expanded case is
		var d = appCtxt.getMsgDialog();
		d.setMessage(ZmMsg.errorSearchNotExpanded);
		d.popup();
		if (expanded) { return; }
	}

	// this method will expand the list depending on the number of email
	// addresses per contact.
	var emailArray = ZmContactsHelper._processSearchResponse(resp, this._includeContactsWithNoEmail); //this._includeContactsWithNoEmail - true in the ZmGroupView case 
	var emailList = AjxVector.fromArray(emailArray);

	if (serverPaginationSupported) {
		this._emailList.addList(emailArray); //this internally calls concat. we do not need "merge" here because we use the _serverContactOffset as a marker of where to search next, never searching a block we already did.
	}
	else {
		this._emailList = emailList;
	}
    var list = this.getSubList();
    if (this.toString() === "ZmContactPicker") {
        list = this.loadSharedGroupContacts(list) || list;
    }
	this._showResults(list);

};

ZmContactPicker.prototype._showResults =
function(aList) {
    var list = aList || this.getSubList();
	// special case 1 - search forward another server block, to fill up a page. Could search several times.
	if (list.size() < ZmContactsApp.SEARCHFOR_MAX && this._serverHasMoreAndPaginationSupported) {
		this.search(null, null, null, null, null, this._serverContactOffset + ZmContactsApp.SEARCHFOR_MAX); //search another page
		return;
	}

	if (this._searchIcon) { //does not exist in ZmGroupView case
		this._searchIcon.className = "";
	}
	this._searchButton.setEnabled(true);

	// special case 2 - no results, and no more to search (that was covered in special case 1) - so display the "no results" text.
	if (list.size() == 0 && this._emailListOffset == 0) {
		this._setResultsInView(list); //empty the list
		this._nextButton.setEnabled(false);
		this._prevButton.setEnabled(false);
		this._setNoResultsHtml();
		return;
	}

	// special case 3 - If the AB ends with a long list of contacts w/o addresses,
	// we may get an empty list.  If that's the case, roll back the offset
	// not 100% sure this case could still happen after all my changes but it was there in the code, so I keep it just in case.
	if (list.size() == 0) {
		this._emailListOffset -= ZmContactsApp.SEARCHFOR_MAX;
		this._emailListOffset  = Math.max(0, this._emailListOffset);
	}

	var more = this._serverHasMoreAndPaginationSupported  //we can get more from the server
				|| (this._emailListOffset + ZmContactsApp.SEARCHFOR_MAX) < this._emailList.size(); //or we have more on the client we didn't yet show
	this._prevButton.setEnabled(this._emailListOffset > 0);
	this._nextButton.setEnabled(more);

	this._resetSearchColHeaders(); // bug #2269 - enable/disable sort column per type of search
	this._setResultsInView(list);
};

ZmContactPicker.prototype.loadSharedGroupContacts =
    function(aList) {

    var listLen,
        listArray,
        contact,
        item,
        i,
        j,
        k,
        sharedContactGroupArray,
        len1,
        jsonObj,
        batchRequest,
        request,
        response;

        listArray = aList.getArray();
        listLen = aList.size();
        sharedContactGroupArray = [];

    for (i = 0 ; i < listLen; i++) {
        item = listArray[i];
        contact = item.__contact;
        if (contact.isGroup() && contact.isShared()) {
            if (this.sharedContactGroups.indexOf(item.value) !== -1) {
                return;
            }
            this.sharedContactGroups.push(item.value);
            sharedContactGroupArray.push(item.value);
        }
    }

    len1 = sharedContactGroupArray.length;
    jsonObj = {BatchRequest:{GetContactsRequest:[],_jsns:"urn:zimbra", onerror:'continue'}};
    batchRequest = jsonObj.BatchRequest;
    request = batchRequest.GetContactsRequest;

    for (j = 0,k =0; j < len1; j++) {
        request.push({ cn: {id: sharedContactGroupArray[j]}, _jsns: 'urn:zimbraMail', derefGroupMember: '1', returnCertInfo: '1', requestId: k++ }); //See ZCS-999 and ZCS-991 for returnCertInfo param
    }
        var respCallback = new AjxCallback(this, this.handleSharedContactResponse,[aList]);
       response =  appCtxt.getAppController().sendRequest({
            jsonObj:jsonObj,
            asyncMode:true,
            callback:respCallback
        });
   };

ZmContactPicker.prototype.handleSharedContactResponse =
    function(aList,response) {

     var contactResponse,
         contactResponseLength,
         listArray,
         listArrayLength,
         sharedGroupMembers,
         i,
         j,
         k,
         resp,
         contact,
         member,
         isGal,
         memberContact,
         loadMember,
         listArrElement,
         sharedGroupMembers;

        if (response && response.getResponse() && response.getResponse().BatchResponse) {
            contactResponse = response.getResponse().BatchResponse.GetContactsResponse;
        }
        if (!contactResponse) {
            return;
        }
        contactResponseLength = contactResponse.length;
        listArray = aList.getArray();
        listArrayLength = aList.size();

        for (k= 0; k < listArrayLength; k++) {
            sharedGroupMembers = [];
            for (j = 0; j < contactResponseLength; j++) {
                resp = contactResponse[j];
                contact = resp.cn[0];

                if (contact.m) {
                    for (i = 0; i < contact.m.length; i++) {
                        member = contact.m[i];
                        isGal = false;
                        if (member.type == ZmContact.GROUP_GAL_REF) {
                            isGal = true;
                        }
                        if (member.cn && member.cn.length > 0) {
                            memberContact = member.cn[0];
                            memberContact.ref = memberContact.ref || (isGal && member.value);
                            loadMember = ZmContact.createFromDom(memberContact, {list: this.list, isGal: isGal});
                            loadMember.isDL = isGal && loadMember.attr[ZmContact.F_type] === "group";
                            appCtxt.cacheSet(member.value, loadMember);
                            listArrElement = listArray[k];
                            if (listArrElement.value === contact.id) {
                                sharedGroupMembers.push( '"'+loadMember.getFullName()+'"' +' <' + loadMember.getEmails() +'>;' ); // Updating the original list with shared members of shared contact group that comes in 'contactResponse'.
                                aList._array[k].address = sharedGroupMembers.join("");
                            }
                        }
                    }
                    ZmContact.prototype._loadFromDom(contact);
                }
            }
        }
        this._showResults(aList); // As async = true, when the response has come, again we render/update the  list with  contactResponse shared contact members,
    };


/**
 * extracted this so it can be used in ZmGroupView where this is different.
 * @param list
 */
ZmContactPicker.prototype._setResultsInView =
function(list) {
	this._chooser.setItems(list);
};

/**
 * extracted this so it can be used in ZmGroupView where this is different.
 * @param list
 */
ZmContactPicker.prototype._setNoResultsHtml =
function(list) {
	this._chooser.sourceListView._setNoResultsHtml();
};


ZmContactPicker.prototype._updateSearchRows =
function(searchFor) {
	var fieldIds = (searchFor == ZmContactsApp.SEARCHFOR_GAL) ? ZmContactPicker.SHOW_ON_GAL : ZmContactPicker.SHOW_ON_NONGAL;
	for (var fieldId in this._searchRow) {
		Dwt.setVisible(this._searchRow[fieldId], AjxUtil.indexOf(fieldIds, fieldId)!=-1);
	}
	for (var fieldId in this._searchField) {
		var field = this._searchField[fieldId];
		if (this._tabGroup.contains(field))
			this._tabGroup.removeMember(field);
	}
	for (var i=0; i<fieldIds.length; i++) {
		this._tabGroup.addMember(this._searchField[fieldIds[i]]);
	}

	this._resizeChooser();
};

/**
 * @private
 */
ZmContactPicker.prototype._handleErrorSearch =
function() {
	this._searchButton.setEnabled(true);
	return false;
};

/**
 * @private
 */
ZmContactPicker.prototype._pageListener =
function(ev) {
	if (ev.item == this._prevButton) {
		this._emailListOffset -= ZmContactsApp.SEARCHFOR_MAX;
		this._emailListOffset  = Math.max(0, this._emailListOffset);
	}
	else {
		this._emailListOffset += ZmContactsApp.SEARCHFOR_MAX;
	}
	this._showResults();
};

/**
 * Gets a sub-list of contacts.
 * 
 * @return	{AjxVector}		a vector of {ZmContact} objects
 */
ZmContactPicker.prototype.getSubList =
function() {
	var size = this._emailList.size();

	var end = this._emailListOffset + ZmContactsApp.SEARCHFOR_MAX;

	if (end > size) {
		end = size;
	}

	var a = (this._emailListOffset < end) ? this._emailList.getArray().slice(this._emailListOffset, end) : [];
	return AjxVector.fromArray(a);
};

/**
 * @private
 */
ZmContactPicker.prototype._searchTypeListener =
function(ev) {
	var oldValue = ev._args.oldValue;
	var newValue = ev._args.newValue;

	if (oldValue != newValue) {
		this._updateSearchRows(newValue);
		this._searchButtonListener();
	}
};

/**
 * @private
 */
ZmContactPicker.prototype._resetSearchColHeaders =
function () {
    var slv = this._chooser.sourceListView;
    var tlv = this._chooser.targetListView;
    slv.headerColCreated = false;
    tlv.headerColCreated = false;
    var isGal = this._searchInSelect && (this._searchInSelect.getValue() == ZmContactsApp.SEARCHFOR_GAL);

    // find the participant column
    var part = 0;
    for (var i = 0; i < slv._headerList.length; i++) {
        var field = slv._headerList[i]._field;
        if (field == ZmItem.F_NAME) {
            part = i;
        }
        if (field == ZmItem.F_DEPARTMENT) {
            slv._headerList[i]._visible = isGal && this._detailedSearch;
        }
    }

    var sortable = isGal ? null : ZmItem.F_NAME;
    slv._headerList[part]._sortable = sortable;
    slv.createHeaderHtml(sortable);

    for (i = 0; i < tlv._headerList.length; i++) {
        if (tlv._headerList[i]._field == ZmItem.F_DEPARTMENT) {
            tlv._headerList[i]._visible = isGal && this._detailedSearch;
        }
    }
    tlv.createHeaderHtml();
};

/**
 * Done choosing addresses, add them to the compose form.
 * 
 * @private
 */
ZmContactPicker.prototype._okButtonListener =
function(ev) {
	var data = this._chooser.getItems();
	DwtDialog.prototype._buttonListener.call(this, ev, [data]);
};

/**
 * Call custom popdown method.
 * 
 * @private
 */
ZmContactPicker.prototype._cancelButtonListener =
function(ev) {
	DwtDialog.prototype._buttonListener.call(this, ev);
	this.popdown();
};

/**
 * @private
 */
ZmContactPicker._keyPressHdlr =
function(ev) {
	var stb = DwtControl.getTargetControl(ev);
	var charCode = DwtKeyEvent.getCharCode(ev);
	if (stb._keyPressCallback && (charCode == 13 || charCode == 3)) {
		stb._keyPressCallback.run();
		return false;
	}
	return true;
};


/***********************************************************************************/

/**
 * Creates a contact chooser.
 * @class
 * This class creates a specialized chooser for the contact picker.
 *
 * @param {DwtComposite}	parent			the contact picker
 * @param {Array}		buttonInfo		transfer button IDs and labels
 * 
 * @extends		DwtChooser
 * 
 * @private
 */
ZmContactChooser = function(params) {
	DwtChooser.call(this, params);
};

ZmContactChooser.prototype = new DwtChooser;
ZmContactChooser.prototype.constructor = ZmContactChooser;

/**
 * @private
 */
ZmContactChooser.prototype._createSourceListView =
function() {
	return new ZmContactChooserSourceListView(this);
};

/**
 * @private
 */
ZmContactChooser.prototype._createTargetListView =
function() {
	return new ZmContactChooserTargetListView(this, (this._buttonInfo.length > 1));
};

/**
 * The item is a AjxEmailAddress. Its address is used for comparison.
 *
 * @param {AjxEmailAddress}	item	an email address
 * @param {AjxVector}	list	list to check in
 * 
 * @private
 */
ZmContactChooser.prototype._isDuplicate =
function(item, list) {
	return list.containsLike(item, item.getAddress);
};

/***********************************************************************************/

/**
 * Creates a source list view.
 * @class
 * This class creates a specialized source list view for the contact chooser.
 * 
 * @param {DwtComposite}	parent			the contact picker
 * 
 * @extends		DwtChooserListView
 * 
 * @private
 */
ZmContactChooserSourceListView = function(parent) {
	DwtChooserListView.call(this, {parent:parent, type:DwtChooserListView.SOURCE});
	this.setScrollStyle(Dwt.CLIP);
};

ZmContactChooserSourceListView.prototype = new DwtChooserListView;
ZmContactChooserSourceListView.prototype.constructor = ZmContactChooserSourceListView;

/**
 * Returns a string representation of the object.
 * 
 * @return		{String}		a string representation of the object
 * @private
 */
ZmContactChooserSourceListView.prototype.toString =
function() {
	return "ZmContactChooserSourceListView";
};

ZmContactChooserSourceListView.prototype.getToolTipContent =
function(ev) {
	
	if (this._hoveredItem) {
		var ttParams = {
			address:		this._hoveredItem.address,
			contact:		this._hoveredItem.__contact,
			ev:				ev,
			noRightClick:	true
		};
		var ttCallback = new AjxCallback(this,
			function(callback) {
				appCtxt.getToolTipMgr().getToolTip(ZmToolTipMgr.PERSON, ttParams, callback);
			});
		return {callback:ttCallback};
	}
	else {
		return "";
	}
};

/**
 * @private
 */
ZmContactChooserSourceListView.prototype._getHeaderList =
function() {
	var headerList = [];
	headerList.push(new DwtListHeaderItem({field:ZmItem.F_TYPE, icon:"Folder", width:ZmMsg.COLUMN_WIDTH_FOLDER_CN}));
	headerList.push(new DwtListHeaderItem({field:ZmItem.F_NAME, text:ZmMsg._name, width:ZmMsg.COLUMN_WIDTH_NAME_CN, resizeable: true}));
	headerList.push(new DwtListHeaderItem({field:ZmItem.F_DEPARTMENT, text:ZmMsg.department, width:ZmMsg.COLUMN_WIDTH_DEPARTMENT_CN, resizeable: true}));
	headerList.push(new DwtListHeaderItem({field:ZmItem.F_EMAIL, text:ZmMsg.email, resizeable: true}));


	return headerList;
};

 
// Override of DwtListView.prototype._resetColWidth to set width; without overrriding causes vertical scrollbars to disapper
// on header resize
ZmContactChooserSourceListView.prototype._resetColWidth =
function() {

	if (!this.headerColCreated) { return; }

	var lastColIdx = this._getLastColumnIndex();
    if (lastColIdx) {
        var lastCol = this._headerList[lastColIdx];
        var lastCell = document.getElementById(lastCol._id);
		if (lastCell) {
			var div = lastCell.firstChild;
			lastCell.style.width = div.style.width = (lastCol._width || ""); 
		}
    }
};

/**
 * override for scrollbars in IE
 * @param headerIdx
 */
ZmContactChooserSourceListView.prototype._calcRelativeWidth =
function(headerIdx) {
	var column = this._headerList[headerIdx];
	if (!column._width || (column._width && column._width == "auto")) {
		var cell = document.getElementById(column._id);
		// UGH: clientWidth is 5px more than HTML-width (20px for IE to deal with scrollbars)
		return (cell) ? (cell.clientWidth - (AjxEnv.isIE ? Dwt.SCROLLBAR_WIDTH : 5)) : null;
	}
	return column._width;
};

/**
 * @private
 */
ZmContactChooserSourceListView.prototype._mouseOverAction =
function(ev, div) {
	DwtChooserListView.prototype._mouseOverAction.call(this, ev, div);
	var id = ev.target.id || div.id;
	var item = this.getItemFromElement(div);
	this._hoveredItem = (id && item) ? item : null;
	return true;
};

/**
 * @private
 */
ZmContactChooserSourceListView.prototype._getCellContents =
function(html, idx, item, field, colIdx, params) {
	if (field == ZmItem.F_EMAIL && AjxEnv.isIE) {
		var maxWidth = AjxStringUtil.getWidth(item.address);
		html[idx++] = "<div style='float; left; overflow: visible; width: " + maxWidth + ";'>";
		idx = ZmContactsHelper._getEmailField(html, idx, item, field, colIdx, params);
		html[idx++] = "</div>";		
	}
	else {
		idx = ZmContactsHelper._getEmailField(html, idx, item, field, colIdx, params);
	}
	return idx;
};

/**
 * Returns a string of any extra attributes to be used for the TD.
 *
 * @param item		[object]	item to render
 * @param field		[constant]	column identifier
 * @param params	[hash]*		hash of optional params
 * 
 * @private
 */
ZmContactChooserSourceListView.prototype._getCellAttrText =
function(item, field, params) {
	if (field == ZmItem.F_EMAIL) {
		return "style='position: relative; overflow: visible;'";
	}
};

/***********************************************************************************/

/**
 * Creates the target list view.
 * @class
 * This class creates a specialized target list view for the contact chooser.
 * 
 * @param {DwtComposite}	parent			the contact picker
 * @param {constant}		showType		the show type
 * @extends		DwtChooserListView
 * 
 * @private
 */
ZmContactChooserTargetListView = function(parent, showType) {
	this._showType = showType; // call before base class since base calls getHeaderList

	DwtChooserListView.call(this, {parent:parent, type:DwtChooserListView.TARGET});

	this.setScrollStyle(Dwt.CLIP);
};

ZmContactChooserTargetListView.prototype = new DwtChooserListView;
ZmContactChooserTargetListView.prototype.constructor = ZmContactChooserTargetListView;

/**
 * Returns a string representation of the object.
 * 
 * @return		{String}		a string representation of the object
 */
ZmContactChooserTargetListView.prototype.toString =
function() {
	return "ZmContactChooserTargetListView";
};

/**
 * @private
 */
ZmContactChooserTargetListView.prototype._getHeaderList =
function() {
	var headerList = [];
	var view = this._view;
	if (this._showType) {
		headerList.push(new DwtListHeaderItem({field:ZmItem.F_TYPE, icon:"ContactsPicker", width:ZmMsg.COLUMN_WIDTH_TYPE_CN}));
	}
	headerList.push(new DwtListHeaderItem({field:ZmItem.F_NAME, text:ZmMsg._name, width:ZmMsg.COLUMN_WIDTH_NAME_CN, resizeable: true}));
    headerList.push(new DwtListHeaderItem({field:ZmItem.F_DEPARTMENT, text:ZmMsg.department, width:ZmMsg.COLUMN_WIDTH_DEPARTMENT_CN, resizeable: true}));
    headerList.push(new DwtListHeaderItem({field:ZmItem.F_EMAIL, text:ZmMsg.email, resizeable: true}));

	return headerList;
};

ZmContactChooserTargetListView.prototype._mouseOverAction =
ZmContactChooserSourceListView.prototype._mouseOverAction;

/**
 * The items are AjxEmailAddress objects.
 * 
 * @private
 */
ZmContactChooserTargetListView.prototype._getCellContents =
function(html, idx, item, field, colIdx, params) {
	if (field == ZmItem.F_TYPE) {
		item.setType(item._buttonId);
		html[idx++] = ZmMsg[item.getTypeAsString()];
		html[idx++] = ":";
	}
	else if (field == ZmItem.F_EMAIL && AjxEnv.isIE) {
		var maxWidth = AjxStringUtil.getWidth(item.address) + 10;
		html[idx++] = "<div style='float; left;  width: " + maxWidth + ";'>";
		idx = ZmContactsHelper._getEmailField(html, idx, item, field, colIdx, params);
		html[idx++] = "</div>";
	}
	else {
		idx = ZmContactsHelper._getEmailField(html, idx, item, field, colIdx);
	}
	return idx;
};


// Override of DwtListView.prototype._resetColWidth to set width; without overrriding causes vertical scrollbars to disapper
// on header resize
ZmContactChooserTargetListView.prototype._resetColWidth =
function() {

	if (!this.headerColCreated) { return; }

	var lastColIdx = this._getLastColumnIndex();

	
    if (lastColIdx) {
        var lastCol = this._headerList[lastColIdx];
        var lastCell = document.getElementById(lastCol._id);
		if (lastCell) {
			var div = lastCell.firstChild;
			lastCell.style.width = div.style.width = (lastCol._width || "");
		}
    }
};

/**
 * override for scrollbars in IE
 * @param headerIdx
 */
ZmContactChooserTargetListView.prototype._calcRelativeWidth =
function(headerIdx) {
	var column = this._headerList[headerIdx];
	if (!column._width || (column._width && column._width == "auto")) {
		var cell = document.getElementById(column._id);
		// UGH: clientWidth is 5px more than HTML-width (20px for IE to deal with scrollbars)
		return (cell) ? (cell.clientWidth - (AjxEnv.isIE ? Dwt.SCROLLBAR_WIDTH : 5)) : null;
	}
	return column._width;
};

/**
 * Returns a string of any extra attributes to be used for the TD.
 *
 * @param item		[object]	item to render
 * @param field		[constant]	column identifier
 * @param params	[hash]*		hash of optional params
 * 
 * @private
 */
ZmContactChooserTargetListView.prototype._getCellAttrText =
function(item, field, params) {
	if (field == ZmItem.F_EMAIL) {
		return "style='position: relative; overflow: visible;'";
	}
};
}
}
if (AjxPackage.define("abook.Contacts")) {
AjxTemplate.register("abook.Contacts#ZmEditContactView", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_DEBUG'></div><table role=\"presentation\" valign='top' width='100%'>";
	buffer[_i++] =  AjxTemplate.expand("#ZmEditContactView_header", data) ;
	buffer[_i++] =  AjxTemplate.expand("#ZmEditContactView_body", data) ;
	buffer[_i++] = "</table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactView"
}, false);
AjxTemplate.register("abook.Contacts", AjxTemplate.getTemplate("abook.Contacts#ZmEditContactView"), AjxTemplate.getParams("abook.Contacts#ZmEditContactView"));

AjxTemplate.register("abook.Contacts#ZmEditContactView_header", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<tr valign='top' class='contactHeaderRow'><td class='contactHeaderCell' align='center'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_IMAGE' tabindex='100'></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_REMOVE_IMAGE_row' style='padding-left:7px;font-size:.8em;white-space:nowrap;'><a id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_VIEW_IMAGE' href='#view' tabindex='101'>";
	buffer[_i++] = ZmMsg.view;
	buffer[_i++] = "</a>\n";
	buffer[_i++] = "\t\t\t\t\t|\n";
	buffer[_i++] = "\t\t\t\t\t<a id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_REMOVE_IMAGE' href='#remove' tabindex='102'>";
	buffer[_i++] = ZmMsg.remove;
	buffer[_i++] = "</a></div></td><td class='contactHeaderCell' valign='bottom'><div style='padding:.125em;'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_FULLNAME'></div></div>";
	 if (appCtxt.get(ZmSetting.TAGGING_ENABLED)) { 
	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_TAG' class='contactHeaderCellRow'></div>";
	 } 
	buffer[_i++] =  AjxTemplate.expand("#ZmEditContactView_headerName", data) ;
	buffer[_i++] =  AjxTemplate.expand("#ZmEditContactView_headerJobInfo", data) ;
	buffer[_i++] = "</td><td class='contactHeaderCell'><table role=\"presentation\" class='ZPropertySheet' cellspacing='6'><tr><td style='text-align:right;white-space:nowrap;'>";
	buffer[_i++] = ZmMsg.fileAsLabel;
	buffer[_i++] = "</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_FILE_AS' tabindex='400'></div></td></tr><tr><td style='text-align:right;white-space:nowrap;'>";
	buffer[_i++] = ZmMsg.locationLabel;
	buffer[_i++] = "</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_FOLDER' tabindex='401'></div></td></tr>";
	 if (appCtxt.multiAccounts) { 
	buffer[_i++] = "<tr><td style='text-align:right;'>";
	buffer[_i++] = ZmMsg.accountLabel;
	buffer[_i++] = "</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_ACCOUNT' tabindex='402'></div></td></tr>";
	 } 
	buffer[_i++] = "</table></td></tr>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactView_header"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactView_headerName", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_NAME_row' class='contactHeaderCellRow'><table role=\"presentation\">";
	 if (appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS)) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_row'><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_PREFIX' class='inlineInput'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_FIRST_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_FIRST' tabindex='190'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_MIDDLE' class='inlineInput'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_MAIDEN' class='inlineInput'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_LAST_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_LAST' tabindex='191'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_SUFFIX' class='inlineInput'></td></tr>";
	 } 
	buffer[_i++] = "<tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PREFIX_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PREFIX' tabindex='200'></div></td><td class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_FIRST' tabindex='201'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_MIDDLE_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_MIDDLE' tabindex='202'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_MAIDEN_row' class='inlineInput'><table role=\"presentation\"><tr><td class='inlineInputWrapper'>(</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_MAIDEN' tabindex='203'></div></td><td class='inlineInputWrapper'>)</td></tr></table></td><td class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_LAST' tabindex='204'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_SUFFIX_row' class='inlineInput'><table role=\"presentation\"><tr><td class='inlineInputWrapper'>,</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_SUFFIX' tabindex='205'></div></td></tr></table></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_DETAILS' tabindex='206'></div></td></tr></table></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_NICKNAME_row' class='contactHeaderCellRow'><table role=\"presentation\"><tr><td class='inlineInputWrapper'>&ldquo;</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_NICKNAME' tabindex='250'></div></td><td class='inlineInputWrapper'>&rdquo;</td></tr></table></div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactView_headerName"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactView_headerJobInfo", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_JOB_row' class='contactHeaderCellRow'><table role=\"presentation\"><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_TITLE_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_TITLE' tabindex='301'></div></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_TITLE_DEPARTMENT_SEP' class='inlineInputWrapper'>&ndash;</td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_DEPARTMENT_row' class='inlineInput'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_DEPARTMENT' tabindex='302'></div></td></tr></table></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_COMPANY_row' class='contactHeaderCellRow'><table role=\"presentation\">";
	 if (appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS)) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_COMPANY_row'><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONETIC_COMPANY' tabindex='304'></div></td></tr>";
	 } 
	buffer[_i++] = "<tr><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_COMPANY' tabindex='305'></div></td></tr></table></div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactView_headerJobInfo"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactView_body", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.emailLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_EMAIL' tabindex='500'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.phoneLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_PHONE' tabindex='600'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.imLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_IM' tabindex='700'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.addressLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_ADDRESS' tabindex='800'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.urlLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_URL' tabindex='900'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.otherLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_OTHER' tabindex='1000'></div></td></tr><tr><td class=rowLabel>";
	buffer[_i++] = ZmMsg.notesLabel;
	buffer[_i++] = "</td><td class=rowValue><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_NOTES' tabindex='1100'></div></td></tr>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactView_body"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewRows", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\"><tbody id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_rows' tabindex='100'></tbody></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewRows"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewDetailsButton", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<span class='ScreenReaderOnly'>";
	buffer[_i++] =  ZmMsg.chooseFields ;
	buffer[_i++] = "</span><span aria-hidden='true' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_title' class='ZButtonBorder ZWidgetTitle'></span>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewDetailsButton"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewInputSelect", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='inputTable'><tr><td class='inputElement'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_input' tabindex='100'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_select' tabindex='200'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewInputSelect"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewInputDoubleSelect", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='inputTable'><tr><td class='inputElement'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_input' tabindex='100'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_select' tabindex='200'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_select2' tabindex='300'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewInputDoubleSelect"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewOther", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='inputTable'><tr><td class='inputElement'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_input' tabindex='100'></div></td><td class='inputElement' style='width:47px;'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_picker' tabindex='200'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_select' tabindex='300'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewOther"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewAddressRow", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='inputTable'><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_row' class='DwtFormRow' valign='top'><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "' tabindex='100'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_add' class='DwtFormRowAdd' tabindex='200'></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_remove' class='DwtFormRowRemove' tabindex='300'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewAddressRow"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewAddressSelect", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='inputTable'><tr valign='top'><td class='inputElement'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_input' tabindex='100' style=\"padding:1px 3px;\"></div></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_select' tabindex='200'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewAddressSelect"
}, false);

AjxTemplate.register("abook.Contacts#ZmEditContactViewAddress", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='addressTable' width='100%'><tr><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_STREET' tabindex='100'></div></td></tr></table><table role=\"presentation\" class='addressTable' width='100%' style='margin-bottom:3px;'><tr><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_CITY' tabindex='200'></div></td><td>,&nbsp;</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_STATE' tabindex='300'></div></td><td>&nbsp;&nbsp;</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_ZIP' tabindex='400'></div></td></tr></table><table role=\"presentation\" class='addressTable' width='100%'><tr><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_COUNTRY' tabindex='500'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmEditContactViewAddress",
	"class": "DwtForm ZmEditContactViewAddress"
}, false);

AjxTemplate.register("abook.Contacts#SplitView", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" width='100%' height='100%' class='ZPropertySheet' cellspacing='6'><col width='230'></col><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_alphabetbar'></td></tr><tr height='100%'><td valign='top'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_splitview' style='position:relative;width:100%;height:100%;'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_listview' style='position:absolute;width:230px;height:100%;left:0;'></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_sash' style='position:absolute;width:8px;height:100%;left:230px;'></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_contentCell' style='position:absolute;left:238px;'><div class='ZmContactInfoView' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_content'></div></div></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_content", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_details' class='contactContentTable'>";
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#SplitView_header", data) ;
	buffer[_i++] =  ZmContactSplitView.showContactEmails(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactPhones(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactIMs(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactAddresses(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactUrls(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactOther(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactNotes(data) ;
	buffer[_i++] =  ZmContactSplitView.showContactDLMembers(data) ;
	buffer[_i++] = "</table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_content"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_group", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	
			var imageUrl = data.imageUrl;
			var defaultImageUrl = data.defaultImageUrl;
			var email = data.email;
			var title = data.title;
			var fullName = data.fullName;
			var phone = data.phone;
			var location = data.location;
			var imgClassName = data.imgClassName || "Person_48";
			var quickAdd = "id='" + data.quickAddId + "'";
			var isEdit = data.isEdit;
			var delId = "id='" + data.delButtonId + "'";
		
	buffer[_i++] = "<table role=\"presentation\" class='contactGroupTable'><tr>";
	 if (imageUrl) { 
	buffer[_i++] = "<td class='contactGroupTableImage' width='48'><img src='";
	buffer[_i++] =  imageUrl ;
	buffer[_i++] = "'\n";
	buffer[_i++] = "\t\t\t";
	 if (AjxEnv.isIE) { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\theight='48'\n";
	buffer[_i++] = "\t\t\t";
	 } else { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\tstyle='max-width:48px;max-height:48px;'\n";
	buffer[_i++] = "\t\t\t";
	 } 
	 if (defaultImageUrl) { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\tonerror='this.onerror=null;this.src=\"";
	buffer[_i++] = data["defaultImageUrl"];
	buffer[_i++] = "\";'\n";
	buffer[_i++] = "\t\t\t";
	 } 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\tborder='0'></td>";
	 } else { 
	buffer[_i++] = "<td class='contactGroupTableImage' width='48'>";
	buffer[_i++] =  AjxImg.getImageHtml(imgClassName) ;
	buffer[_i++] = "</td>";
	 } 
	buffer[_i++] = "<td class='contactGroupTableContent'>";
	 if (fullName) {
	buffer[_i++] = "<div><strong>";
	buffer[_i++] = fullName;
	buffer[_i++] = "</strong></div>";
	 } 
	 if (title) { 
	buffer[_i++] = "<div><span>";
	buffer[_i++] = title;
	buffer[_i++] = "</span></div>";
	}
	 if (email) {
	buffer[_i++] = "<div>";
	buffer[_i++] = email;
	 if (data.isInline) {
	buffer[_i++] = AjxImg.getImageHtml("Plus", "", quickAdd);
	 }
	buffer[_i++] = "</div>";
	}
	 if (phone) {
	buffer[_i++] = "<div>";
	buffer[_i++] = phone;
	buffer[_i++] = "</div>";
	}
	 if (location) {
	buffer[_i++] = "<div>";
	buffer[_i++] = location;
	buffer[_i++] = "</div>";
	}
	buffer[_i++] = "</td>";
	 if (isEdit) {
	buffer[_i++] = "<td width='36' valign='top'>";
	buffer[_i++] = AjxImg.getImageHtml("Delete", "", delId);
	buffer[_i++] = "</td>";
	}
	buffer[_i++] = "</tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_group"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_header", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	
			var contact = data.contact;
			var imageUrl = contact.getImageUrl();
			var imageLabel = contact.isDistributionList() ? ZmMsg.distributionList : contact.isGroup() ? ZmMsg.group : ZmMsg.contact;
			var defaultImageUrl = data.defaultImageUrl;
			var accountName = appCtxt.multiAccounts && contact.account && contact.account.getDisplayName();
			var tagsId = data.id + '_tags_contact'; 
			var imgClassName = contact.isDistributionList() ? 'Group_48' : contact.isGroup() ? 'GroupPerson_48' : 'Person_48';
			var isEdit = data.isEdit;
		
	buffer[_i++] = "<tr class='headerRow' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_headerRow'>";
	 if (imageUrl) { 
	buffer[_i++] = "<td class='rowLabel'><div class='headerRowImage'><img src='";
	buffer[_i++] =  imageUrl ;
	buffer[_i++] = "'\n";
	buffer[_i++] = "\t\t\t\t";
	 if (AjxEnv.isIE) { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\t\theight='48'\n";
	buffer[_i++] = "\t\t\t\t";
	 } else { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\t\tstyle='max-width:48px;max-height:48px;'\n";
	buffer[_i++] = "\t\t\t\t";
	 } 
	 if (defaultImageUrl) { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\t\tonerror='this.onerror=null;this.src=\"";
	buffer[_i++] = data["defaultImageUrl"];
	buffer[_i++] = "\";'\n";
	buffer[_i++] = "\t\t\t\t";
	 } 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\t\tborder=0 alt='";
	buffer[_i++] =  imageLabel ;
	buffer[_i++] = "'></div></td>";
	 } else { 
	buffer[_i++] = "<td class='rowLabel'><div class='headerRowImage'>";
	buffer[_i++] =  AjxImg.getImageHtml({
							imageName: imgClassName, altText: imageLabel
						}) ;
	buffer[_i++] = "</div></td>";
	 } 
	buffer[_i++] = "<td class='rowValue'>";
	buffer[_i++] =  AjxTemplate.expand("#SplitView_headerName", data) ;
	buffer[_i++] =  AjxTemplate.expand("#SplitView_headerJobInfo", data) ;
	 if (data.dlInfo && data.dlInfo.description) { 
	buffer[_i++] = "<div>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.dlInfo.description) ;
	buffer[_i++] = "</div>";
	 } 
	 if (!isEdit) {
	buffer[_i++] = "<div class='contactHeaderTags' id='";
	buffer[_i++] =  tagsId ;
	buffer[_i++] = "'></div>";
	 } 
	buffer[_i++] = "</td><td class='rowType' valign='top'>";
	 if (accountName) { 
	buffer[_i++] = "<table role=\"presentation\" class='ZPropertySheet' cellspacing='6'><tr><td class='contactLocation'>";
	buffer[_i++] =  ZmMsg.accountLabel ;
	buffer[_i++] = "</td><td class='companyFolder'>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(accountName) ;
	buffer[_i++] = "</td></tr></table>";
	 } 
	 if (!isEdit) { 
	 if (data.addrbook) { 
	buffer[_i++] =  AjxImg.getImageHtml({
								imageName: data.addrbook.getIconWithColor(),
								altText: ZmMsg.locationLabel
							}) ;
	buffer[_i++] = "&nbsp;";
	buffer[_i++] =  data.addrbook.getName() ;
	 } else { 
	buffer[_i++] =  AjxImg.getImageHtml('GAL') ;
	buffer[_i++] = "&nbsp;";
	buffer[_i++] =  ZmMsg.GAL ;
	 } 
	 } else { 
	buffer[_i++] = "<table role=\"presentation\" class='ZPropertySheet' cellspacing='6'><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_locationLabel' class='contactLocation'>";
	buffer[_i++] =  ZmMsg.locationLabel ;
	buffer[_i++] = "</td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_LOCATION_FOLDER'></td></tr></table>";
	 } 
	buffer[_i++] = "</td></tr>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_header"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_headerName", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	
			var contact = data.contact;
			var isPhonetic = appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS);
			var fullnameHtml = contact.getFullNameForDisplay(isPhonetic);
			if (!isPhonetic) {
				fullnameHtml = AjxStringUtil.htmlEncode(fullnameHtml);
			}
			var nickname = contact.getAttr(ZmContact.F_nickname);
			var isEdit = data.isEdit;
			if (isEdit && data.contact.isDistributionList()) {
				fullnameHtml = '';
			}
		
	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_contactName' class='contactHeader ";
	buffer[_i++] =  data.isInTrash ? "Trash" : "" ;
	buffer[_i++] = "'>";
	buffer[_i++] =  fullnameHtml ;
	buffer[_i++] = "</div>";
	 if (nickname && !isEdit) { 
	buffer[_i++] = "<div class='contactHeader'>&ldquo;";
	buffer[_i++] =  AjxStringUtil.htmlEncode(nickname) ;
	buffer[_i++] = "&rdquo;</div>";
	 } 
	 else if (isEdit) {
	buffer[_i++] = "<div class='contactHeader2'>";
		if (data.contact.isDistributionList()) {
				
	buffer[_i++] = "<table role=\"presentation\"><tr><td class='contactHeader2'>";
	
						if (data.usernameEditable) { 
	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_groupNameParent'></div>";
			} else { 
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.username);
			}
				
	buffer[_i++] = "</td><td class='contactHeader2'>&nbsp;@\n";
	buffer[_i++] = "\t\t\t";
	
						if (data.domainEditable) { 
	buffer[_i++] = "<input type='text' autocomplete='off' size='20' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_groupNameDomain' class='Text'>";
			} else { 
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.domain);
			} 
	buffer[_i++] = "</td></tr></table>";
	 } else { 
	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_groupNameParent'></div>";
	 } 
	buffer[_i++] = "</div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_tags' style='padding-top:5px;'></div>";
	 } 

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_headerName"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_headerJobInfo", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	
			var contact = data.contact;
			var isPhonetic = appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS);
			var company = contact.getAttr(ZmContact.F_company);
			var title = contact.getAttr(ZmContact.F_jobTitle);
			var department = contact.getAttr(ZmContact.F_department);
		
	 if (title || department) { 
	buffer[_i++] = "<div class='contactHeader2'>";
	 if (title) { 
	buffer[_i++] = "<span>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(title) ;
	buffer[_i++] = "</span>";
	 } 
	 if (department) { 
	 if (title) { 
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t\t&nbsp;&ndash;&nbsp;\n";
	buffer[_i++] = "\t\t\t";
	 } 
	buffer[_i++] = "<span>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(department) ;
	buffer[_i++] = "</span>";
	 } 
	buffer[_i++] = "</div>";
	 } 
	 if (company) {
			var phoneticCompany = isPhonetic && contact.getAttr(ZmContact.F_phoneticCompany);
			var companyHtml = AjxStringUtil.htmlRubyEncode(company, phoneticCompany);
			
	buffer[_i++] = "<div class='contactHeader2'>";
	buffer[_i++] =  companyHtml;
	buffer[_i++] = "</div>";
	 } 

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_headerJobInfo"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_dlmember-collapsed", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<td ";
	buffer[_i++] =  data.expandTdText ;
	buffer[_i++] = ">";
	buffer[_i++] =  data.html.join("") ;
	buffer[_i++] = "</td><td class=rowLabel>";
	buffer[_i++] =  data.contact.canExpand ? ZmMsg.membersLabel : "" ;
	buffer[_i++] = "</td><td colspan='2' class='rowValue'>";
	buffer[_i++] =  data.value ;
	buffer[_i++] = "</td>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_dlmember-collapsed"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_dlmember-expanded", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<td ";
	buffer[_i++] =  data.first ? data.expandTdText : ""  ;
	buffer[_i++] = ">";
	buffer[_i++] =  data.html.join("") ;
	buffer[_i++] = "</td><td class='rowLabel'>";
	buffer[_i++] =  data.first ? ZmMsg.membersLabel : "" ;
	buffer[_i++] = "</td><td colspan='2' class='rowValue'>";
	buffer[_i++] =  data.value ;
	buffer[_i++] = "</td>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_dlmember-expanded"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_address_value", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<tr valign='top'><td class='rowLabel'><!-- only show the first label visually --><label for='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "' class='";
	buffer[_i++] =  data.first ? "" : "ScreenReaderOnly" ;
	buffer[_i++] = "'>\n";
	buffer[_i++] = "\t\t\t\t";
	buffer[_i++] = data["label"];
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t</label></td><td id='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "' class='rowValue'>";
	 var address = data.address;
				if (address.Street) { 
	buffer[_i++] = "<div>";
	buffer[_i++] = address.Street;
	buffer[_i++] = "</div>";
	 }
				if (address.City || address.State || address.PostalCode) { 
	buffer[_i++] = "<div><span>";
	buffer[_i++] = address.City;
	buffer[_i++] = "</span>";
	 if (address.City && address.State) { 
	buffer[_i++] = ", ";
	 } 
	buffer[_i++] = "<span>";
	buffer[_i++] = address.State;
	buffer[_i++] = "</span>";
	 if (address.City || address.State) { 
	buffer[_i++] = "&nbsp;";
	 } 
	buffer[_i++] = "<span>";
	buffer[_i++] = address.PostalCode;
	buffer[_i++] = "</span></div>";
	 }
				if (address.Country) { 
	buffer[_i++] = "<div>";
	buffer[_i++] = address.Country;
	buffer[_i++] = "</div>";
	 } 
	buffer[_i++] = "</td><td class='rowType'><label for='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "'>\n";
	buffer[_i++] = "\t\t\t\t&nbsp;";
	buffer[_i++] =  data.type ;
	buffer[_i++] = "</label></td></tr>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_address_value"
}, false);

AjxTemplate.register("abook.Contacts#SplitView_list_item", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<tr valign='top'>";
	 if (data.isDL) { 
	buffer[_i++] = "<td><div class='ImgBlank_16'></div></td>";
	 } 
	buffer[_i++] = "<td class='rowLabel' id='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_label'><!-- only show the first label visually --><label for='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_label' class='";
	buffer[_i++] =  data.first ? "" : "ScreenReaderOnly" ;
	buffer[_i++] = "'>\n";
	buffer[_i++] = "\t\t\t";
	buffer[_i++] = data["label"];
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t</label></td><td id='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_value' ";
	buffer[_i++] =  data.type ? "" : "colspan='2'" ;
	buffer[_i++] = " class='rowValue' aria-labelledby='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_label ";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_type'>";
	buffer[_i++] =  data.value ;
	buffer[_i++] = "</td>";
	 if (data.type) { 
	buffer[_i++] = "<td id='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_type' class='rowType'><label for='";
	buffer[_i++] = data["labelId"];
	buffer[_i++] = "_label'>\n";
	buffer[_i++] = "\t\t\t\t&nbsp;";
	buffer[_i++] = data["type"];
	buffer[_i++] = "\n";
	buffer[_i++] = "\t\t\t</label></td>";
	 } 
	buffer[_i++] = "</tr>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitView_list_item"
}, false);

AjxTemplate.register("abook.Contacts#SplitViewGroup", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table class='contactContentTable'>";
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#SplitView_header", data) ;
	buffer[_i++] = "</table>";
	 if (data.dlInfo) { 
	buffer[_i++] = "<table role=\"presentation\" class='contactHeaderSubTable NoneBg ZPropertySheet' width='100%' cellspacing='6'><tr><td><ul id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_subscriptionMsg' style='margin:.5em 0;padding-left:53px;'></ul></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_subscriptionButton'></td></tr></table>";
	 } 
	buffer[_i++] = "<div class='contactGroupList' style='overflow:auto;'>";
	buffer[_i++] =  ZmContactSplitView.showContactGroup(data) ;
	buffer[_i++] = "</div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SplitViewGroup"
}, false);

AjxTemplate.register("abook.Contacts#SimpleView-NoResults", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" width='100%'><tr><td class='NoResults'>";
	buffer[_i++] =  AjxMsg.noResults ;
	buffer[_i++] = "</td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#SimpleView-NoResults"
}, false);

AjxTemplate.register("abook.Contacts#GroupView", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table class='ZPropertySheet'>";
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#SplitView_header", data) ;
	buffer[_i++] = "</table>";
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#GroupViewMembers", data) ;

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#GroupView"
}, false);

AjxTemplate.register("abook.Contacts#DlView", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table class='ZPropertySheet'>";
	buffer[_i++] =  AjxTemplate.expand('abook.Contacts#SplitView_header', data) ;
	buffer[_i++] = "</table><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_tabViewContainer'></div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#DlView"
}, false);

AjxTemplate.register("abook.Contacts#DlPropertiesView", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	 if (data.contact.isDistributionList()) { 
	buffer[_i++] = "<table role=\"presentation\" width='100%' class='ZmEditContactView'><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.displayNameLabel ;
	buffer[_i++] = "</td><td class='rowValue'><div class='DwtInputField'><input type='text' autocomplete='off' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlDisplayName' class='Text' style='width:45em;'></div></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.descriptionLabel ;
	buffer[_i++] = "</td><td class='rowValue'><div class='DwtInputField'><textarea rows='4' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlDesc' style='width:45em;'></textarea></div></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.directoryLabel ;
	buffer[_i++] = "</td><td class='rowValue rowCheckGroup'><input type='checkbox' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlHideInGal'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlHideInGal'>";
	buffer[_i++] =  ZmMsg.hiddenInGlobal ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.mailLabel ;
	buffer[_i++] = "</td><td class='rowValue rowRadioGroup'><input type='radio' name='mailPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyANYONE'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyANYONE'>";
	buffer[_i++] =  ZmMsg.mailAnyone ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='mailPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyMEMBERS'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyMEMBERS'>";
	buffer[_i++] =  ZmMsg.mailOnlyMembers ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='mailPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyINTERNAL'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicyINTERNAL'>";
	buffer[_i++] =  ZmMsg.mailOnlyInternal ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='mailPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicySPECIFIC'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlMailPolicySPECIFIC'>";
	buffer[_i++] =  ZmMsg.mailOnlySpecific ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><div class='DwtInputField'><input type='text' autocomplete='off' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlListSpecificMailers' class='Text' style='margin-left:2em;width:43em;'></div></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.newSubscriptionsLabel ;
	buffer[_i++] = "</td><td class='rowValue rowRadioGroup'><input type='radio' name='subscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyACCEPT'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyACCEPT'>";
	buffer[_i++] =  ZmMsg.policyAccept ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='subscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyAPPROVAL'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyAPPROVAL'>";
	buffer[_i++] =  ZmMsg.policyApproval ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='subscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyREJECT'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlSubscriptionPolicyREJECT'>";
	buffer[_i++] =  ZmMsg.policyReject ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.unsubscriptionsRequestLabel ;
	buffer[_i++] = "</td><td class='rowValue rowRadioGroup' ><input type='radio' name='unsubscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyACCEPT'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyACCEPT'>";
	buffer[_i++] =  ZmMsg.policyAccept ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='unsubscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyAPPROVAL'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyAPPROVAL'>";
	buffer[_i++] =  ZmMsg.policyApproval ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'></td><td class='rowValue'><input type='radio' name='unsubscriptionPolicy' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyREJECT'><label for='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlUnsubscriptionPolicyREJECT'>";
	buffer[_i++] =  ZmMsg.policyReject ;
	buffer[_i++] = "</label></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.listOwnersLabel ;
	buffer[_i++] = "</td><td class='rowValue'><div class='DwtInputField'><input type='text' autocomplete='off' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlListOwners' class='Text' style='width:45em;'></div></td></tr><tr><td class='rowLabel'>";
	buffer[_i++] =  ZmMsg.notesLabel ;
	buffer[_i++] = "</td><td class='rowValue'><div class='DwtInputField'><textarea rows='4' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_dlNotes' style='width:45em;'></textarea></div></td></tr></table>";
	 } 

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#DlPropertiesView"
}, false);

AjxTemplate.register("abook.Contacts#GroupViewMembers", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" width='100%' height='100%' class='ZPropertySheet' cellspacing='6'><tr><!-- content: left pane --><td width='50%' valign='top'><div class='groupMembers' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_groupMembers'></div></td><!-- content: right pane --><td width='50%' valign='top'><table role=\"presentation\" width='100%'><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchFieldsRow'><td width='50'>&nbsp;</td><td><table role=\"presentation\" class='ZPropertySheet' cellspacing='6'>";
	 if (data.detailed) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameRow'><td align='right'>";
	buffer[_i++] =  ZmMsg.nameLabel ;
	buffer[_i++] = "</td><td><input type='text' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameField'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchButton' valign='top' rowspan='3'></td></tr>";
	 if (appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS)) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticRow'><td align='right'>";
	buffer[_i++] =  ZmMsg.phoneticNameLabel ;
	buffer[_i++] = "</td><td><input type='text' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticField'></td></tr>";
	 } 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailRow'><td align='right'>";
	buffer[_i++] =  ZmMsg.emailAddrLabel ;
	buffer[_i++] = "</td><td><input type='text' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailField'></td></tr><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentRow'><td align='right'>";
	buffer[_i++] =  ZmMsg.departmentLabel ;
	buffer[_i++] = "</td><td><input type='text' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentField'></td></tr>";
	 } else { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchRow'><td align='right'>";
	buffer[_i++] =  ZmMsg.findLabel ;
	buffer[_i++] = "</td><td><input type='text' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchField'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchButton'></td></tr>";
	 } 
	 if (data.showSearchIn) { 
	buffer[_i++] = "<tr><td align='right'>";
	buffer[_i++] =  ZmMsg.searchIn ;
	buffer[_i++] = "</td><td colspan='2' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_listSelect'></td></tr>";
	 } 
	buffer[_i++] = "</table></td></tr><tr><td width='50'><table role=\"presentation\" width='100%' class='ZPropertySheet' cellspacing='6'><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_addButton'></td></tr><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_addAllButton'></td></tr></table></td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_listView'></div><table role=\"presentation\" width='100%' class='ZPropertySheet' cellspacing='6'><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_navButtonsRow'><td width='100%'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_prevButton'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_nextButton'></td></tr></table></td></tr><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_manualAddRow'><td width='50'><table role=\"presentation\" width='100%' class='ZPropertySheet' cellspacing='6'><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_addNewButton'></td></tr></table></td><td><table role=\"presentation\" width='100%' class='ZPropertySheet' cellspacing='6'><tr><td>";
	buffer[_i++] = ZmMsg.enterAddrBelow;
	buffer[_i++] = "</td></tr><tr><td><textarea class='groupMembers' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_addNewField'></textarea></td></tr></table></td></tr></table></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#GroupViewMembers"
}, false);

AjxTemplate.register("abook.Contacts#ZmAlphabetBar", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<center><table role=\"presentation\" id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_alphabet' width='80%' class='AlphabetBarTable'' cellspacing='0'><tr>";
	 for (var i = 0; i < data.numLetters; i++) { 
	buffer[_i++] = "<td _idx=\"";
	buffer[_i++] =  i ;
	buffer[_i++] = "\"\trole='button'\n";
	buffer[_i++] = "\t\t\t\t\t\tclass='DwtButton AlphabetBarCell'";
	 if (data.alphabet[i].length > 1) { 
	buffer[_i++] = " colspan='";
	buffer[_i++] =  data.alphabet[i].length/2 + 1 ;
	buffer[_i++] = "' ";
	 } 
	buffer[_i++] = " >";
	buffer[_i++] =  data.alphabet[i] ;
	buffer[_i++] = "</td>";
	 } 
	buffer[_i++] = "</tr></table></center>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmAlphabetBar"
}, false);

AjxTemplate.register("abook.Contacts#Tooltip", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" width='250'><tr><td colspan='2' valign='top'><div style='border-bottom:1px solid black;margin-bottom:2px;'><table role=\"presentation\" width='100%'><tr valign='bottom'><td style='font-weight:bold;'>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.entryTitle); ;
	buffer[_i++] = "</td><td align='right'>";
	buffer[_i++] =  AjxImg.getImageHtml(data.contact.getIcon()); ;
	buffer[_i++] = "</td></tr></table></div></td></tr>";
	 if (data.contact.isGroup()) {
				var members = data.contact.getGroupMembers().good.getArray();
				for (var i = 0; i < members.length; i++) {
			
	buffer[_i++] = "<tr><td width=20>";
	buffer[_i++] =  AjxImg.getImageHtml("Message") ;
	buffer[_i++] = "</td><td>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(members[i].toString()) ;
	buffer[_i++] = "</td></tr>";
		} 
			} else {
			
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"fullName"}) ;
	 var fields = ZmMsg.contactTooltipWorkInfoOrder.split(",");
				   var i;
				   for (i = 0; i < fields.length; i++) {
				
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:fields[i]}) ;
	
				   }
				
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"mobilePhone"}) ;
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"workPhone"}) ;
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"homePhone"}) ;
	 if (data.contact.isGal) {
					var emails = data.contact.getEmails();
					for (var i = 0; i < emails.length && emails.length <= 3; i++) {
				
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"email", val:emails[i]}) ;
		} 
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"notes", truncate:500}) ;
		} else { 
	buffer[_i++] =  AjxTemplate.expand("abook.Contacts#AddTooltipEntry", {data:data, field:"email", val:data.contact._lookupEmail || data.email}) ;
		} 
		} 
	 if (data.hint) { 
	buffer[_i++] = "<tr><td colspan='2'><div class='TooltipHint'><hr color=black size=1>";
	buffer[_i++] =  data.hint ;
	buffer[_i++] = "</div></td></tr>";
	 } 
	buffer[_i++] = "</table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#Tooltip"
}, false);

AjxTemplate.register("abook.Contacts#TooltipNotInAddrBook", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	 if (data.hint) { 
	buffer[_i++] = "<table role=\"presentation\" width='250'><tr><td valign='top'>";
	 if (appCtxt.get(ZmSetting.CONTACTS_ENABLED)) { 
	buffer[_i++] = "<div style='border-bottom:1px solid black;'><table role=\"presentation\" width='100%'><tr valign='bottom'><td style='font-weight:bold;'>";
	buffer[_i++] =  ZmMsg.newContact ;
	buffer[_i++] = "</td><td align='right'>";
	buffer[_i++] =  AjxImg.getImageHtml("NewContact"); ;
	buffer[_i++] = "</td></tr></table></div>";
	 } 
	buffer[_i++] = "<div class='TooltipNotInAddrBook'>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.addrstr) ;
	buffer[_i++] = "</div><hr color='black' size='1'><div class='TooltipHint'>";
	buffer[_i++] =  data.hint ;
	buffer[_i++] = "</div></td></tr></table>";
	 } else { 
	buffer[_i++] = "<span style='font-weight:bold;'> ";
	buffer[_i++] =  ZmMsg.emailLabel ;
	buffer[_i++] = " </span>";
	buffer[_i++] =  AjxStringUtil.htmlEncode(data.addrstr) ;
	 } 

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#TooltipNotInAddrBook"
}, false);

AjxTemplate.register("abook.Contacts#AddTooltipEntry", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	 if (data.val == null) {
			var isPhonetic = appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS);
			data.val = data.field == "fullName" ? data.data.contact.getFullNameForDisplay(isPhonetic) : data.data.contact.getAttr(data.field);
		} 
	 if (data.val != null && data.val != "") { 
	buffer[_i++] = "<tr valign='top'><td class='ZmTooltipLabel'>";
	buffer[_i++] =  AjxMessageFormat.format(ZmMsg.makeLabel,
											 AjxStringUtil.htmlEncode(ZmContact._AB_FIELD[data.field])) ;
	buffer[_i++] = "</td>";
	 if (data.truncate) { 
	buffer[_i++] = "<td>";
	buffer[_i++] =  AjxMessageFormat.format(ZmMsg.makeLabel,
											AjxStringUtil.htmlEncode(data.val.substring(0,data.truncate))) ;
	buffer[_i++] = "</td>";
	 } else {
				var field = "";
				if (data.field == "fullName" ) {
					field = isPhonetic ? data.val : AjxStringUtil.htmlEncode(data.val); 
				}
				else {
					field = AjxStringUtil.htmlEncode(data.val);
				}
			
	buffer[_i++] = "<td style='white-space:nowrap;'>";
	buffer[_i++] =  field ;
	buffer[_i++] = "</td>";
	 } 
	buffer[_i++] = "</tr>";
	 } 

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#AddTooltipEntry"
}, false);

AjxTemplate.register("abook.Contacts#ZmContactPicker", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<div class='ZmContactPicker'><table role=\"presentation\" width='100%'><tr><td width='600'><table role=\"presentation\" class='ZPropertySheet' cellspacing='6' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchTable'>";
	 if (data.detailed) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameRow'><td align='right'><span id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameFieldLbl'>";
	buffer[_i++] =  ZmMsg.nameLabel ;
	buffer[_i++] = "</span></td><td><input type='text' placeholder=\"";
	buffer[_i++] =  ZmMsg.contactPickerHint ;
	buffer[_i++] = "\" autocomplete='off' size='30' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameField' aria-labelledby='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchNameFieldLbl'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchButton' rowspan='3' valign='top'></td><td width='20' valign='middle'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchIcon' class='ImgSearch'></div></td></tr>";
	 if (appCtxt.get(ZmSetting.PHONETIC_CONTACT_FIELDS)) { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticRow'><td align='right'><span id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticNameFieldLbl'>";
	buffer[_i++] =  ZmMsg.phoneticNameLabel ;
	buffer[_i++] = "</span></td><td><input type='text' placeholder=\"";
	buffer[_i++] =  ZmMsg.contactPickerPhoneticHint ;
	buffer[_i++] = "\" autocomplete='off' size='30' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticField' aria-labelledby='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchPhoneticNameFieldLbl'></td></tr>";
	 } 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailRow'><td align='right'><span id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailFieldLbl'>";
	buffer[_i++] =  ZmMsg.emailAddrLabel ;
	buffer[_i++] = "</span></td><td><input type='text' placeholder=\"";
	buffer[_i++] =  ZmMsg.contactPickerEmailHint ;
	buffer[_i++] = "\"autocomplete='off' size='30' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailField' aria-labelledby='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchEmailFieldLbl'></td></tr><tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentRow'><td align='right'><span id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentFieldLbl'>";
	buffer[_i++] =  ZmMsg.departmentLabel ;
	buffer[_i++] = "</span></td><td><input type='text' placeholder=\"";
	buffer[_i++] =  ZmMsg.contactPickerDepartmentHint ;
	buffer[_i++] = "\" autocomplete='off' size='30' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentField' aria-labelledby='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchDepartmentFieldLbl'></td></tr>";
	 } else { 
	buffer[_i++] = "<tr id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchRow'><td><input type='text' placeholder=\"";
	buffer[_i++] =  ZmMsg.contactPickerHint ;
	buffer[_i++] = "\" autocomplete='off' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchField' style='width:25em;' aria-label='";
	buffer[_i++] = ZmMsg.search;
	buffer[_i++] = "'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchButton'></td><td width='20' valign='middle'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchIcon' class='ImgSearch'></div></td></tr>";
	 } 
	buffer[_i++] = "</table></td><td align='right' valign='bottom'>";
	 if (data.showSelect) { 
	buffer[_i++] = "<table role=\"presentation\" class='ZPropertySheet' cellspacing='6'><tr><td class='Label nobreak' id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_listSelectLbl'>";
	buffer[_i++] =  ZmMsg.showNames ;
	buffer[_i++] = "</td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_listSelect'></td></tr></table>";
	 } 
	buffer[_i++] = "</td></tr></table><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_chooser'></div><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_paging' style='margin-left:5px;'><table role=\"presentation\" width='100%'><tr><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_pageLeft'></td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_pageRight' align='right'></td></tr></table></div></div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmContactPicker"
}, false);

AjxTemplate.register("abook.Contacts#ZmContactSearch", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_preamble'></div><table role=\"presentation\" width='100%'><tr><td><table role=\"presentation\"><tr><td width='20' valign='middle'><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchIcon' class='ImgSearch'></div></td><td><input type='text' autocomplete='off' size='30' nowrap id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchField'>&nbsp;</td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_searchButton'></td></tr></table></td><td align='right'>";
	 if (data.showSelect) { 
	buffer[_i++] = "<table role=\"presentation\"><tr><td class='Label nobreak'>";
	buffer[_i++] =  ZmMsg.showNames ;
	buffer[_i++] = "&nbsp;</td><td id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_folders'></td></tr></table>";
	 } 
	buffer[_i++] = "</td></tr></table><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_results' class='contactSearchResultsDiv'></div>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#ZmContactSearch"
}, false);

AjxTemplate.register("abook.Contacts#QuickAddPrompt", 
function(name, params, data, buffer) {
	var _hasBuffer = Boolean(buffer);
	data = (typeof data == "string" ? { id: data } : data) || {};
	buffer = buffer || [];
	var _i = buffer.length;

	buffer[_i++] = "<table role=\"presentation\" class='ZPropertySheet' cellspacing='6'><tr><td class='ZmFieldLabelRight'>";
	buffer[_i++] = ZmMsg.firstNameLabel;
	buffer[_i++] = "</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_FIRST_NAME'></div></td></tr><tr><td class='ZmFieldLabelRight'>";
	buffer[_i++] = ZmMsg.lastNameLabel;
	buffer[_i++] = " </td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_LAST_NAME'></div></td></tr><tr><td class='ZmFieldLabelRight'>";
	buffer[_i++] = ZmMsg.emailLabel;
	buffer[_i++] = "</td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_EMAIL'></div></td></tr><tr><td class='ZmFieldLabelRight'>";
	buffer[_i++] = ZmMsg.addressBookLabel;
	buffer[_i++] = " </td><td><div id='";
	buffer[_i++] = data["id"];
	buffer[_i++] = "_ADDR_BOOK'></div></td></tr></table>";

	return _hasBuffer ? buffer.length : buffer.join("");
},
{
	"id": "abook.Contacts#QuickAddPrompt"
}, false);

}
